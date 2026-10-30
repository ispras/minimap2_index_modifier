// Copyright 2023 karpulevich
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

#include <htslib/vcf.h>
#include "mmpriv.h"
#include "minimap.h"
#include "linked_vcf_list.h"
#include "rmin.h"

struct node *head = NULL;
struct node *current = NULL;

static int gt_error(bcf_hdr_t *hdr, bcf1_t *rec, int sample, const char *reason)
{
    const char *chrom = rec ? bcf_hdr_id2name(hdr, rec->rid) : "header";
    fprintf(stderr, "[ERROR::GT] %s:%ld sample=%s: %s\n",
            chrom ? chrom : "unknown",
            rec ? (long)rec->pos + 1 : 0L,
            sample >= 0 ? hdr->samples[sample] : "-", reason);
    return -1;
}

/* Ploidy is fixed for each sample, while HTSlib's padded row width may vary. */
static int decode_gt(bcf_hdr_t *hdr, bcf1_t *rec, int *ploidy, char *copies)
{
    int32_t *gt = NULL;
    int capacity = 0, n, offset = 0, status = -1;
    int samples = bcf_hdr_nsamples(hdr);
    if (samples <= 0) return gt_error(hdr, rec, -1, "GT requires samples");
    bcf_unpack(rec, BCF_UN_STR);
    if (rec->n_allele != 2) return gt_error(hdr, rec, -1, "only biallelic records are supported");
    for (int i = 0; i < 2; ++i)
        if (!*rec->d.allele[i] || strspn(rec->d.allele[i], "ACGTNacgtn") != strlen(rec->d.allele[i]))
            return gt_error(hdr, rec, -1, "unsupported non-sequence allele");
    n = bcf_get_genotypes(hdr, rec, &gt, &capacity);
    if (n <= 0 || n % samples != 0) {
        gt_error(hdr, rec, -1, "missing or malformed FORMAT/GT");
        goto done;
    }
    for (int i = 0; i < samples; ++i) {
        int width = n / samples, count = 0;
        int32_t *row = gt + (size_t)i * width;
        while (count < width && row[count] != bcf_int32_vector_end) ++count;
        if (count < 1 || count > 2) {
            gt_error(hdr, rec, i, "only haploid or diploid GT is supported"); goto done;
        }
        for (int j = count; j < width; ++j)
            if (row[j] != bcf_int32_vector_end) {
                gt_error(hdr, rec, i, "invalid GT vector padding"); goto done;
            }
        for (int j = 0; j < count; ++j) {
            if (bcf_gt_is_missing(row[j])) {
                gt_error(hdr, rec, i, "missing allele is not a REF allele"); goto done;
            }
            if (bcf_gt_allele(row[j]) < 0 || bcf_gt_allele(row[j]) > 1) {
                gt_error(hdr, rec, i, "unsupported allele index"); goto done;
            }
        }
        if (count == 2 && bcf_gt_allele(row[0]) != bcf_gt_allele(row[1]) && !bcf_gt_is_phased(row[1])) {
            gt_error(hdr, rec, i, "unphased heterozygous GT has no defined carrier phase"); goto done;
        }
        if (ploidy[i] && ploidy[i] != count) {
            gt_error(hdr, rec, i, "sample ploidy changes between records"); goto done;
        }
        ploidy[i] = count;
        for (int j = 0; j < count; ++j) {
            if (copies) copies[offset] = '0' + bcf_gt_allele(row[j]);
            ++offset;
        }
    }
    status = offset;
done:
    free(gt);
    return status;
}

/* Return actual sample copies before deduplication, including real REF copies. */
int mm_gt_decode_window(bcf_hdr_t *hdr, struct node *first, struct node *end,
                        char **carriers, int *rows, int *stride)
{
    int samples = bcf_hdr_nsamples(hdr), sites = 0, count = 0;
    int *ploidy = NULL;
    char *copies = NULL, *matrix = NULL;
    *carriers = NULL; *rows = *stride = 0;
    if (samples <= 0 || samples > INT_MAX / 2)
        return gt_error(hdr, first ? first->rec : NULL, -1, "GT requires a bounded nonempty sample set (samples)");
    for (struct node *v = first; v != end; v = v->next) {
        if (!v || sites == INT_MAX - 1)
            return gt_error(hdr, NULL, -1, "invalid or oversized variant window");
        ++sites;
    }
    if (!sites) return 0;
    if ((size_t)(sites + 1) > SIZE_MAX / ((size_t)samples * 2) ||
        (size_t)sites > SIZE_MAX / sizeof(unsigned long) || (size_t)sites > SIZE_MAX / sizeof(char *))
        return gt_error(hdr, NULL, -1, "carrier matrix size overflow");
    ploidy = (int *)calloc(samples, sizeof(*ploidy));
    copies = (char *)malloc((size_t)samples * 2);
    matrix = (char *)calloc((size_t)samples * 2, (size_t)sites + 1);
    if (!ploidy || !copies || !matrix) {
        gt_error(hdr, NULL, -1, "out of memory decoding GT"); goto fail;
    }
    int site = 0;
    for (struct node *v = first; v != end; v = v->next, ++site) {
        count = decode_gt(hdr, v->rec, ploidy, copies);
        if (count < 0) goto fail;
        for (int i = 0; i < count; ++i) matrix[(size_t)i * (sites + 1) + site] = copies[i];
    }
    free(ploidy); free(copies);
    *carriers = matrix; *rows = count; *stride = sites + 1;
    return 0;
fail:
    free(ploidy); free(copies); free(matrix);
    return -1;
}

/* Validate the complete input before any index part can be written. */
int mm_gt_validate_file(const char *filename)
{
    htsFile *fp = bcf_open(filename, "r");
    bcf_hdr_t *hdr = NULL;
    bcf1_t *rec = NULL;
    int *ploidy = NULL, result = -1, ret;
    if (!fp || !(hdr = bcf_hdr_read(fp))) {
        fprintf(stderr, "[ERROR::GT] cannot read VCF header: %s\n", filename); goto done;
    }
    int samples = bcf_hdr_nsamples(hdr);
    if (samples <= 0 || samples > INT_MAX / 2) {
        gt_error(hdr, NULL, -1, "GT requires a bounded nonempty sample set"); goto done;
    }
    ploidy = (int *)calloc(samples, sizeof(*ploidy));
    rec = bcf_init();
    if (!ploidy || !rec) { gt_error(hdr, NULL, -1, "out of memory validating GT"); goto done; }
    while ((ret = bcf_read(fp, hdr, rec)) == 0) {
        if (rec->errcode) { gt_error(hdr, rec, -1, "malformed VCF record"); goto done; }
        if (decode_gt(hdr, rec, ploidy, NULL) < 0) goto done;
    }
    if (ret != -1) { gt_error(hdr, rec, -1, "malformed VCF record"); goto done; }
    result = 0;
done:
    free(ploidy);
    if (rec) bcf_destroy(rec);
    if (hdr) bcf_hdr_destroy(hdr);
    if (fp && bcf_close(fp) != 0) result = -1;
    return result;
}

// display the list
void printList(){
    struct node *p = head;
    printf("\n[");

    //start from the beginning
    while(p != NULL) {
        printf(" %lu ",p->pos);
        p = p->next;
    }
    printf("]");
}

// display the GT list
void printGTList(bcf_hdr_t *hdr){
    struct node *p = head;


    //start from the beginning
    while(p != NULL) {
        int32_t *gt_arr = NULL, ngt_arr = 0;
        int ngt;
        ngt = bcf_get_genotypes(hdr, p->rec, &gt_arr, &ngt_arr);
        if ( ngt > 0 && bcf_hdr_nsamples(hdr) > 0 ) {
            int i, nsmpl = bcf_hdr_nsamples(hdr);
            printf("%s %lu REF:%s ALT:%s\n", bcf_hdr_id2name(hdr, p->CHR_ID), p->pos, p->REF, p->ALT);

            int max_ploidy = ngt/nsmpl;
            for (i=0; i<nsmpl; i++)
            {
                int32_t *ptr = gt_arr + i*max_ploidy;
                for (int j = 0; j < max_ploidy && ptr[j] != bcf_int32_vector_end; ++j) {
                    if (j) putchar(bcf_gt_is_phased(ptr[j]) ? '|' : '/');
                    if (bcf_gt_is_missing(ptr[j])) putchar('.');
                    else printf("%d", bcf_gt_allele(ptr[j]));
                }
                putchar('\t');
            }
            printf("\n");
        }
        free(gt_arr);

        p = p->next;
    }
}

int ifexists(char* z, const int u, const int s, char* v)
{
    for (int i = 0; i < u; i++)
        if (!strncmp(&z[(size_t)i * s], v, s)) {
            return 1;
        }
    return 0;
}

int calculate_haplotypes(mm_idx_t * mi, bcf_hdr_t *hdr, struct node *window_start_pointer, struct node *current_pointer, unsigned long curr_pos, mm128_v *p) {
    // bcf1_t::pos is 0-based; VCF POS is pos + 1
    struct node *c_pointer = current_pointer;
    struct node *w_start_pointer = window_start_pointer;

    int MAX_SNP = 10;

    char *gt_array = NULL;
    int arr_i, arr_j, arr_len;

    if (!(mi->flag & MM_PARSE_HT)) {
        int snp_num = 0;
        while (w_start_pointer != c_pointer) {
            w_start_pointer = w_start_pointer->next;
            snp_num += 1;
	    if (snp_num > MAX_SNP) return 0;
        }
        /* no variant in the window: nothing to do, and 1 << (snp_num - 1) below would be UB */
        if (snp_num == 0) return 0;
        arr_i = (int) pow(2, snp_num);
        arr_j = snp_num + 1;
        gt_array = (char *) calloc(arr_i * arr_j, sizeof(char));
        if (!gt_array) return gt_error(hdr, NULL, -1, "out of memory decoding combinations");
        for (unsigned i = 0; i < (unsigned) pow(2, snp_num); i++) {
            int k = 0;
            for (int j = 1 << (snp_num - 1); j > 0; j = j / 2) {
                gt_array[i * arr_j + k] = (i & j) ? '1' : '0';
                k++;
            }
        }
        arr_len = arr_i;
    }
    else {
        if (mm_gt_decode_window(hdr, w_start_pointer, c_pointer, &gt_array, &arr_i, &arr_j) < 0)
            return -1;
        if (!arr_i) return 0;
        MAX_SNP = arr_j - 1;

        //remove non unique strings
        int k = 1;
        for (int i = 1; i < arr_i; i++) {
            if (!ifexists(gt_array, k, arr_j, &gt_array[(size_t)i * arr_j])) {
                if (i != k) memcpy(&gt_array[(size_t)k * arr_j], &gt_array[(size_t)i * arr_j], arr_j * sizeof(char));
                k++;
            }
        }
        arr_len = k;
    }

    //prepare SNP combinations
    //Array format:
    //REF_arr char array - REF (for control)
    //ALT_arr char array - ALT
    //POS_all ulong array - positions
    //CHR - chromosome
    //N_SNP - length
    const char * CHR = bcf_hdr_id2name(hdr, window_start_pointer->CHR_ID);

    char **REF_arr = (char **)malloc((size_t)MAX_SNP * sizeof(*REF_arr));
    char **ALT_arr = (char **)malloc((size_t)MAX_SNP * sizeof(*ALT_arr));
    unsigned long *POS_all = (unsigned long *)malloc((size_t)MAX_SNP * sizeof(*POS_all));
    if (!gt_array || !REF_arr || !ALT_arr || !POS_all) {
        fprintf(stderr, "[ERROR::GT] out of memory emitting variant carriers\n");
        free(gt_array); free(REF_arr); free(ALT_arr); free(POS_all);
        return -1;
    }
    rmin_rec_begin(CHR, curr_pos + 1);
    for(int i = 0; i<arr_len; i++) {
        struct node *local_c_pointer = current_pointer;
        struct node *local_w_start_pointer = window_start_pointer;
        int N_SNP = 0;
        int local_snp_num = 0;

        while (local_w_start_pointer != local_c_pointer) {
            if(gt_array[(size_t)i * arr_j + local_snp_num] == '1') {
                REF_arr[N_SNP] = local_w_start_pointer->REF;
                ALT_arr[N_SNP] = local_w_start_pointer->ALT;
                POS_all[N_SNP] = (unsigned long)(local_w_start_pointer->pos + 1);
                N_SNP += 1;
            }

            local_snp_num += 1;
            local_w_start_pointer = local_w_start_pointer->next;
        }
        if (N_SNP > 0) {
            add_variants(mi, CHR, REF_arr, ALT_arr, POS_all, N_SNP, curr_pos + 1, p);
        }
    }
    rmin_rec_end();
    free(REF_arr); free(ALT_arr); free(POS_all);
    free(gt_array);
    return 0;
}

int handleGTList(mm_idx_t * mi, bcf_hdr_t *hdr, mm128_v *p){
    struct node *window_start_pointer = head;
    struct node *current_pointer = head;
    struct node *window_end_pointer = head;

    int gap = mi->k - 1;

    if (head == NULL)
        return 0;

    while(window_end_pointer->next != NULL)
    {
        while(window_end_pointer != NULL && (current_pointer->pos < (unsigned long)gap || window_end_pointer->pos > (current_pointer->pos - gap))){
            window_end_pointer = window_end_pointer->next;
        }

        while(window_start_pointer->pos > (current_pointer->pos + gap)){
            window_start_pointer = window_start_pointer->next;
        }

        //calculate start and end pointers
        if (calculate_haplotypes(mi, hdr, window_start_pointer, window_end_pointer, current_pointer->pos, p) < 0) return -1;

        current_pointer = current_pointer->next;
        window_end_pointer = current_pointer;
    }

    //last SNP batch
    while(window_start_pointer->pos > current_pointer->pos + gap){
        window_start_pointer = window_start_pointer->next;
    }

    return calculate_haplotypes(mi, hdr, window_start_pointer, window_end_pointer->next, current_pointer->pos, p);
}

//insertion at the beginning
int insertatbegin(unsigned long pos, bcf1_t *rec, int CHR_ID, char * REF, char * ALT){
    //create a link
    struct node *lk = (struct node*) malloc(sizeof(struct node));
    if (!lk) {
        fprintf(stderr, "[ERROR::GT] out of memory storing variant record\n");
        return -1;
    }
    lk->pos = pos;
    lk->rec = rec;
    lk->CHR_ID = CHR_ID;
    lk->REF = REF;
    lk->ALT = ALT;

    // point it to old first node
    lk->next = head;

    //point first to new first node
    head = lk;
    return 0;
}

void insertatend(unsigned long pos, bcf1_t *rec, int CHR_ID, char * REF, char * ALT){
    //create a link
    struct node *lk = (struct node*) malloc(sizeof(struct node));
    lk->pos = pos;
    lk->rec = rec;
    lk->CHR_ID = CHR_ID;
    lk->REF = REF;
    lk->ALT = ALT;

    struct node *linkedlist = head;

    // point it to old first node
    while(linkedlist->next != NULL)
        linkedlist = linkedlist->next;

    //point first to new first node
    linkedlist->next = lk;
}
void insertafternode(struct node *list, unsigned long pos, bcf1_t *rec, int CHR_ID, char * REF, char * ALT){
    struct node *lk = (struct node*) malloc(sizeof(struct node));
    lk->pos = pos;
    lk->rec = rec;
    lk->CHR_ID = CHR_ID;
    lk->REF = REF;
    lk->ALT = ALT;

    lk->next = list->next;
    list->next = lk;
}

void deleteatbegin(){
    struct node *temp = head;
    head = temp->next;
    free(temp->ALT);
    free(temp->REF);
    bcf_destroy(temp->rec);
    free(temp);
}

void deleteList(){
    while (head != NULL)
    {
        deleteatbegin();
    }
}

bool isListEmpty(){
    return head == NULL;
}

void deletenode(int key){
    struct node *temp = head, *prev;
    if (temp != NULL && temp->pos == key) {
        head = temp->next;
        free(temp);
        return;
    }

    // Find the key to be deleted
    while (temp != NULL && temp->pos != key) {
        prev = temp;
        temp = temp->next;
    }

    // If the key is not present
    if (temp == NULL) return;

    // Remove the node
    prev->next = temp->next;
    free(temp);
}

int searchlist(int key){
    struct node *temp = head;
    while(temp != NULL) {
        if (temp->pos == key) {
            return 1;
        }
        temp=temp->next;
    }
    return 0;
}
