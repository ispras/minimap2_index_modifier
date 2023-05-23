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
#include <stdbool.h>

#include <htslib/vcf.h>
#include "mmpriv.h"
#include "minimap.h"

struct node {
   unsigned long pos;
   bcf1_t *rec;
   int CHR_ID;
   char * REF;
   char * ALT;
   struct node *next;
};
struct node *head = NULL;
struct node *current = NULL;

// display the list
void printList(){
   struct node *p = head;
   printf("\n[");

   //start from the beginning
   while(p != NULL) {
      printf(" %d ",p->pos);
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
      if ( ngt > 0 ) {
         int i, j, nsmpl = bcf_hdr_nsamples(hdr);
         printf("%s %d REF:%s ALT:%s\n", bcf_hdr_id2name(hdr, p->CHR_ID), p->pos, p->REF, p->ALT);

         int max_ploidy = ngt/nsmpl;
         for (i=0; i<nsmpl; i++)
         {
               int32_t *ptr = gt_arr + i*max_ploidy;
               printf("%d|%d\t",bcf_gt_allele(ptr[0]), bcf_gt_allele(ptr[1]));
         }
         printf("\n");
         free(gt_arr);
      }

      p = p->next;
   }
}

int ifexists(char* z[], int u, char* v)
{
    int i;
    for (i=0; i<u;i++)
        if (!strcmp(z[i],v)) return (1);
    return (0);
}


void calculate_haplotypes(mm_idx_t * mi, bcf_hdr_t *hdr, struct node *window_start_pointer, struct node *current_pointer){
   //Even though vcf uses 1-based indexing (i.e. first base is base 1), htslib internally uses 0-based indexing (i.e. bcf1_t::pos is 0 based).
   //http://wresch.github.io/2014/11/18/process-vcf-file-with-htslib.html
   //so we need use (pos + 1)
   struct node *c_pointer = current_pointer;
   struct node *w_start_pointer = window_start_pointer;


   int MAX_SNP = 40;
   int MAX_HAPLOTYPES = 10000;
   char gt_array[MAX_HAPLOTYPES][MAX_SNP];
   memset(gt_array,'\0', MAX_HAPLOTYPES*MAX_SNP*sizeof(char));
   int sample_num = bcf_hdr_nsamples(hdr);
   int snp_num = 0;


   while (w_start_pointer != c_pointer) {
      int32_t *gt_arr = NULL, ngt_arr = 0;
      int ngt = bcf_get_genotypes(hdr, w_start_pointer->rec, &gt_arr, &ngt_arr);

      if ( ngt > 0 ) {
         int max_ploidy = ngt/sample_num;

         for (int i=0; i<sample_num; i++)
         {
            int32_t *ptr = gt_arr + i*max_ploidy;
            gt_array[i*max_ploidy][snp_num] = bcf_gt_allele(ptr[0]) + '0';
            gt_array[i*max_ploidy+1][snp_num] = bcf_gt_allele(ptr[1]) + '0';
         }
         free(gt_arr);
      }

      //printf("%d ", w_start_pointer->pos+1);
      w_start_pointer = w_start_pointer->next;
      snp_num += 1;
   }

   //printf("\n");


   char * tmp_array [MAX_HAPLOTYPES];

   //remove non unique strings
   //https://subscription.packtpub.com/book/programming/9781838641108/1/ch01lvl1sec06/finding-the-unique-elements-in-an-array
   tmp_array[0]=gt_array[0];
   int k=1;
   int i;
   for (i = 1;i < sample_num*2;i++)
   {
      if(!ifexists(tmp_array,k,gt_array[i]))
      {
         tmp_array[k]=gt_array[i];
         k++;
      }
   }

   //prepare SNP combinations
   //Array format:
   //REF_arr char array - REF (for control)
   //ALT_arr char array - ALT
   //POS_all ulong array - positions
   //CHR - chromosome
   //N_SNP - length
   char * CHR = bcf_hdr_id2name(hdr, window_start_pointer->CHR_ID);

   for(int i = 0; i<k; i++) {
      //printf("%s\n", tmp_array[i]);
      struct node *local_c_pointer = current_pointer;
      struct node *local_w_start_pointer = window_start_pointer;
      int N_SNP = 0;
      int local_snp_num = 0;

      char * REF_arr[MAX_SNP];
      char * ALT_arr[MAX_SNP];
      unsigned long * POS_all[MAX_SNP];

      while (local_w_start_pointer != local_c_pointer) {


         if(tmp_array[i][local_snp_num] == '1') {
            //printf("CHR=%s ALT=%s POS=%d ", bcf_hdr_id2name(hdr, local_w_start_pointer->CHR_ID), local_w_start_pointer->ALT, local_w_start_pointer->pos + 1);
            REF_arr[N_SNP] = local_w_start_pointer->REF;
            ALT_arr[N_SNP] = local_w_start_pointer->ALT;
            POS_all[N_SNP] = (unsigned long)(local_w_start_pointer->pos + 1);
            N_SNP += 1;
         }

         local_snp_num += 1;
         local_w_start_pointer = local_w_start_pointer->next;
      }
      if (N_SNP > 0) {
         add_variants(mi, CHR, REF_arr, ALT_arr, POS_all, N_SNP);
         //printf("-\n");
      }
   }

   //printf("==============\n");

   //print all SNP combinations
   // for(int i = 0; i<10000; i++) {
   //    if(strcmp(gt_array[i], "\0")) {
   //       printf("%s\n", gt_array[i]);
   //    }
   // }

   //print all SNP combinations
   // for(int i = 0; i<10000; i++) {
   //    if(strcmp(gt_array[i], "\0")) {
   //       printf("%s\n", gt_array[i]);
   //    }
   // }


   // for (int j = 0; j < sample_num; j++){
   //    printf("%s\n", gt_array[j*2]);
   //    printf("%s\n", gt_array[j*2+1]);
   // }

   // for (int i = 0; i < snp_num; i++){
   //    for (int j = 0; j < sample_num; j++){
   //       //printf("%c|", gt_array[j*2][i] );
   //       //printf("%c \t", gt_array[j*2+1][i] );
   //    }
   //    //printf("\n");
   // }


   //printf("CHR = %s;\n\n", bcf_hdr_id2name(hdr, window_start_pointer->CHR_ID));
   return;
}

void hadleGTList(mm_idx_t * mi, bcf_hdr_t *hdr, int window, int window_shift){
   struct node *current_pointer = head;
   struct node *window_start_pointer = head;

   if (head == NULL)
      return;

   int window_start_pos = head->pos;
   int window_chr_id_start = head->CHR_ID;

   while(current_pointer->next != NULL)
   {
      current_pointer = current_pointer->next;

      if(window_chr_id_start != current_pointer->CHR_ID) {
         calculate_haplotypes(mi, hdr, window_start_pointer, current_pointer);
         window_start_pointer = current_pointer;
         window_chr_id_start = current_pointer->CHR_ID;
         window_start_pos = current_pointer->pos;
      }

      //calculate and shift
      if((int)(current_pointer->pos) <= (window_start_pos - window)) {
         calculate_haplotypes(mi, hdr, window_start_pointer, current_pointer);

         int new_start_pos = window_start_pos - window_shift;
         if(new_start_pos <=0) break;

         //window shift
         while (window_start_pointer->pos > new_start_pos){
            window_start_pointer = window_start_pointer->next;
         }

         window_chr_id_start = window_start_pointer->CHR_ID;
         window_start_pos = window_start_pointer->pos;
         current_pointer = window_start_pointer;
      }
   }
   //last SNP batch
   calculate_haplotypes(mi, hdr, window_start_pointer, current_pointer->next);
}

//insertion at the beginning
void insertatbegin(unsigned long pos, bcf1_t *rec, int CHR_ID, char * REF, char * ALT){

   //create a link
   struct node *lk = (struct node*) malloc(sizeof(struct node));
   lk->pos = pos;
   lk->rec = rec;
   lk->CHR_ID = CHR_ID;
   lk->REF = REF;
   lk->ALT = ALT;

   // point it to old first node
   lk->next = head;

   //point first to new first node
   head = lk;
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
   printf("DELETE_LIST\n");
   while (head != NULL)
   {
      deleteatbegin();
   }
   printf("DELETE_LIST2\n");
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