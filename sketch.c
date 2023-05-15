#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#define __STDC_LIMIT_MACROS
#include "kvec.h"
#include "mmpriv.h"

unsigned char seq_nt4_table[256] = {
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

static inline uint64_t hash64(uint64_t key, uint64_t mask)
{
	key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
	key = key ^ key >> 28;
	key = (key + (key << 31)) & mask;
	return key;
}

typedef struct { // a simplified version of kdq
	int front, count;
	int a[32];
} tiny_queue_t;

static inline void tq_push(tiny_queue_t *q, int x)
{
	q->a[((q->count++) + q->front) & 0x1f] = x;
}

static inline int tq_shift(tiny_queue_t *q)
{
	int x;
	if (q->count == 0) return -1;
	x = q->a[q->front++];
	q->front &= 0x1f;
	--q->count;
	return x;
}

/**
 * Find symmetric (w,k)-minimizers on a DNA sequence
 *
 * @param km     thread-local memory pool; using NULL falls back to malloc()
 * @param str    DNA sequence
 * @param len    length of $str
 * @param w      find a minimizer for every $w consecutive k-mers
 * @param k      k-mer size
 * @param rid    reference ID; will be copied to the output $p array
 * @param is_hpc homopolymer-compressed or not
 * @param p      minimizers
 *               p->a[i].x = kMer<<8 | kmerSpan
 *               p->a[i].y = rid<<32 | lastPos<<1 | strand
 *               where lastPos is the position of the last base of the i-th minimizer,
 *               and strand indicates whether the minimizer comes from the top or the bottom strand.
 *               Callers may want to set "p->n = 0"; otherwise results are appended to p
 */
void mm_sketch(void *km, const char *str, int len, int w, int k, uint32_t rid, int is_hpc, mm128_v *p)
{
	uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, kmer[2] = {0,0};
	int i, j, l, buf_pos, min_pos, kmer_span = 0;
	mm128_t buf[256], min = { UINT64_MAX, UINT64_MAX };
	tiny_queue_t tq;

	assert(len > 0 && (w > 0 && w < 256) && (k > 0 && k <= 28)); // 56 bits for k-mer; could use long k-mers, but 28 enough in practice
	memset(buf, 0xff, w * 16);
	memset(&tq, 0, sizeof(tiny_queue_t));
	kv_resize(mm128_t, km, *p, p->n + len/w);

	for (i = l = buf_pos = min_pos = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)str[i]];
		mm128_t info = { UINT64_MAX, UINT64_MAX };
		if (c < 4) { // not an ambiguous base
			int z;
			if (is_hpc) {
				int skip_len = 1;
				if (i + 1 < len && seq_nt4_table[(uint8_t)str[i + 1]] == c) {
					for (skip_len = 2; i + skip_len < len; ++skip_len)
						if (seq_nt4_table[(uint8_t)str[i + skip_len]] != c)
							break;
					i += skip_len - 1; // put $i at the end of the current homopolymer run
				}
				tq_push(&tq, skip_len);
				kmer_span += skip_len;
				if (tq.count > k) kmer_span -= tq_shift(&tq);
			} else kmer_span = l + 1 < k? l + 1 : k;
			kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
			kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
			if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
			z = kmer[0] < kmer[1]? 0 : 1; // strand
			++l;
			if (l >= k && kmer_span < 256) {
				info.x = hash64(kmer[z], mask) << 8 | kmer_span;
				info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | z;
			}
		} else l = 0, tq.count = tq.front = 0, kmer_span = 0;
		buf[buf_pos] = info; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below
		if (l == w + k - 1 && min.x != UINT64_MAX) { // special case for the first window - because identical k-mers are not stored yet
			for (j = buf_pos + 1; j < w; ++j)
				if (min.x == buf[j].x && buf[j].y != min.y) kv_push(mm128_t, km, *p, buf[j]);
			for (j = 0; j < buf_pos; ++j)
				if (min.x == buf[j].x && buf[j].y != min.y) kv_push(mm128_t, km, *p, buf[j]);
		}
		if (info.x <= min.x) { // a new minimum; then write the old min
			if (l >= w + k && min.x != UINT64_MAX) kv_push(mm128_t, km, *p, min);
			min = info, min_pos = buf_pos;
		} else if (buf_pos == min_pos) { // old min has moved outside the window
			if (l >= w + k - 1 && min.x != UINT64_MAX) kv_push(mm128_t, km, *p, min);
			for (j = buf_pos + 1, min.x = UINT64_MAX; j < w; ++j) // the two loops are necessary when there are identical k-mers
				if (min.x >= buf[j].x) min = buf[j], min_pos = j; // >= is important s.t. min is always the closest k-mer
			for (j = 0; j <= buf_pos; ++j)
				if (min.x >= buf[j].x) min = buf[j], min_pos = j;
			if (l >= w + k - 1 && min.x != UINT64_MAX) { // write identical k-mers
				for (j = buf_pos + 1; j < w; ++j) // these two loops make sure the output is sorted
					if (min.x == buf[j].x && min.y != buf[j].y) kv_push(mm128_t, km, *p, buf[j]);
				for (j = 0; j <= buf_pos; ++j)
					if (min.x == buf[j].x && min.y != buf[j].y) kv_push(mm128_t, km, *p, buf[j]);
			}
		}
		if (++buf_pos == w) buf_pos = 0;
	}
	if (min.x != UINT64_MAX)
		kv_push(mm128_t, km, *p, min);
}


void mm_idx_manipulate(mm_idx_t * mi, char *vcf_with_variants) {
    FILE *vcf;
    vcf = fopen(vcf_with_variants, "r");

    if(vcf == NULL) {
        perror("Error in opening VCF file");
        exit(-1);
    }

    int VCF_READ_BUFFER_SIZE = 128;
    char buf[VCF_READ_BUFFER_SIZE];
    while (fgets(buf, VCF_READ_BUFFER_SIZE - 2, vcf)) {
        char *token = strtok(buf, "\t");
        char *str[5];
        // Read line
        for (int i = 0; i < 5; i++) {
            str[i] = token;
            token = strtok(NULL, "\t");
        }

        if (str[0][0] == '\n' || str[0][0] == '\r') return;

        const char *snp_contig_name = str[0];
        const char *snp_position = str[1];
        // Should be used for testing
        //const char * snp_from = str[2];
        // WARNING! str[4][0] because of \n as str[4][1]
        const char snp_to = str[4][0];

        //Find seq
        uint64_t contig_offset;
        uint64_t seq_num;
        for (int i = 0; i < mi->n_seq; i++) {
            if (strcmp(snp_contig_name, mi->seq[i].name) == 0) {
                contig_offset = mi->seq[i].offset;
                seq_num = i;
            }
        }

        int SIDE_SIZE = (mi->k - 1) + mi->w;
        // Calculate number of chunks:
        // side chunks: take k-mer size, subtract 1 and add window size
        // divided by chunk size and multiplied by 2 as it has 2 sides, and one for center
        int SEQ_CHUNK_NUMBER = SIDE_SIZE / 8 * 2 + 1;
        // add extra two side chunks if (mi->k - 1 + 10) is not a multiple of 8
        int EXTRA_GAP = (8 - SIDE_SIZE % 8) % 8;
        SEQ_CHUNK_NUMBER = (EXTRA_GAP) ? SEQ_CHUNK_NUMBER + 2 : SEQ_CHUNK_NUMBER;
        uint32_t seq[SEQ_CHUNK_NUMBER];

        for (int i = 0; i < SEQ_CHUNK_NUMBER; i++) {
            if (
                    // Out of bounds
                    (contig_offset == 0 && (atol(snp_position) - 1) / 8 - (SEQ_CHUNK_NUMBER / 2) + i < 0) ||
                    ((contig_offset + atol(snp_position) - 1) / 8 - (SEQ_CHUNK_NUMBER / 2) + i) * 8 >= contig_offset + mi->seq[seq_num].len ||
                    ((contig_offset + atol(snp_position) - 1) / 8 - (SEQ_CHUNK_NUMBER / 2) + i + 1) * 8 <= contig_offset
                    )
                seq[i] = 1145324612; // ALL N
            else {
                seq[i] = mi->S[(contig_offset + atol(snp_position) - 1) / 8 - (SEQ_CHUNK_NUMBER / 2) + i];
                // At left bound
                if (((contig_offset + atol(snp_position) - 1) / 8 - (SEQ_CHUNK_NUMBER / 2) + i) * 8 < contig_offset &&
                    ((contig_offset + atol(snp_position) - 1) / 8 - (SEQ_CHUNK_NUMBER / 2) + i + 1) * 8 > contig_offset) {
                    if (contig_offset % 8 == 0)
                        seq[i] = 1145324612; // ALL N
                    else {
                        seq[i] = seq[i] >> (4 * (contig_offset % 8));
                        for (int j = 0; j < contig_offset % 8; j++)
                            seq[i] = (seq[i] << 4) + 4;
                    }
                }
                // At right bound
                if (((contig_offset + atol(snp_position) - 1) / 8 - (SEQ_CHUNK_NUMBER / 2) + i) * 8 < contig_offset + mi->seq[seq_num].len &&
                    ((contig_offset + atol(snp_position) - 1) / 8 - (SEQ_CHUNK_NUMBER / 2) + i + 1) * 8 > contig_offset + mi->seq[seq_num].len) {
                    seq[i] = seq[i] << (4 * (8 - (contig_offset + mi->seq[seq_num].len) % 8));
                    for (int j = 0; j < 8 - (contig_offset + mi->seq[seq_num].len) % 8; j++)
                        seq[i] = (seq[i] >> 4) | 1073741824; // FIRST N
                }
            }
        }

        char original_ref_seq[SEQ_CHUNK_NUMBER * 8 + 1];
        original_ref_seq[SEQ_CHUNK_NUMBER * 8] = '\0';
        for (int i = 0; i < SEQ_CHUNK_NUMBER; i++) {
            uint32_t tmp_seq = seq[i];
            for (int j = 0; j < 8; j++) {
                uint32_t tmp = tmp_seq % 16;
                switch (tmp) {
                    case 0:
                        original_ref_seq[i * 8 + j] = 'A';
                        break;
                    case 1:
                        original_ref_seq[i * 8 + j] = 'C';
                        break;
                    case 2:
                        original_ref_seq[i * 8 + j] = 'G';
                        break;
                    case 3:
                        original_ref_seq[i * 8 + j] = 'T';
                        break;
                    case 4:
                        original_ref_seq[i * 8 + j] = 'N';
                }
                tmp_seq = tmp_seq / 16;
            }
        }
        // str[3] - REF
        // str[4] - ALT
        if (strlen(str[3]) == 1 && strlen(str[4]) == 1) {
            //SNP
            char new_ref_seq[SEQ_CHUNK_NUMBER * 8 + 1];
            memcpy(new_ref_seq, original_ref_seq, SEQ_CHUNK_NUMBER * 8 + 1);
            new_ref_seq[EXTRA_GAP + SIDE_SIZE + (contig_offset + atol(snp_position) - 1) % 8] = snp_to;

            //Finds minimizer in window
            mm128_v minimizer_array = {0, 0, 0};
            mm_sketch(0, &new_ref_seq[EXTRA_GAP + (contig_offset + atol(snp_position) - 1) % 8], SIDE_SIZE * 2 + 1, mi->w, mi->k,
                      0, mi->flag & MM_I_HPC, &minimizer_array);

            for (int i = 0; i < minimizer_array.n; i++) {
                if (minimizer_array.a[i].y < SIDE_SIZE * 2) continue;
                if (minimizer_array.a[i].y > (SIDE_SIZE + mi->k) * 2 - 1) continue;
                minimizer_array.a[i].y = (seq_num << 32) + (atol(snp_position) - SIDE_SIZE - 1 + minimizer_array.a[i].y / 2) * 2 +
                                         (minimizer_array.a[i].y % 2);
                mm_idx_push(mi, minimizer_array.a[i].x, minimizer_array.a[i].y);
            }
        } else if (strlen(str[3]) > 1 && strlen(str[4]) == 1) {
            //Deletion
            //Gather extra seq
            int ext_chunk_count = (strlen(str[3]) - 2) / 8 + 1;
            char * original_ref_seq_ext;
            original_ref_seq_ext = malloc(sizeof(char) * (8 * ext_chunk_count + 1));
            original_ref_seq_ext[8 * ext_chunk_count] = '\0';
            for (int i = 0; i < ext_chunk_count; i++) {
                uint32_t tmp_seq;
                for (int i = 0; i < ext_chunk_count; i++) {
                    if (
                        // Out of bounds
                            (contig_offset == 0 && (atol(snp_position) - 1) / 8 + (SEQ_CHUNK_NUMBER / 2) + i + 1 < 0) ||
                            ((contig_offset + atol(snp_position) - 1) / 8 + (SEQ_CHUNK_NUMBER / 2) + i + 1) * 8 >= contig_offset + mi->seq[seq_num].len ||
                            ((contig_offset + atol(snp_position) - 1) / 8 + (SEQ_CHUNK_NUMBER / 2) + i + 2) * 8 <= contig_offset
                            )
                        tmp_seq = 1145324612; // ALL N
                    else {
                        tmp_seq = mi->S[(contig_offset + atol(snp_position) - 1) / 8 + (SEQ_CHUNK_NUMBER / 2) + i + 1];
                        // At left bound
                        if (((contig_offset + atol(snp_position) - 1) / 8 + (SEQ_CHUNK_NUMBER / 2) + i + 1) * 8 < contig_offset &&
                            ((contig_offset + atol(snp_position) - 1) / 8 + (SEQ_CHUNK_NUMBER / 2) + i + 2) * 8 > contig_offset) {
                            if (contig_offset % 8 == 0)
                                tmp_seq = 1145324612; // ALL N
                            else {
                                tmp_seq = tmp_seq >> (4 * (contig_offset % 8));
                                for (int j = 0; j < contig_offset % 8; j++)
                                    tmp_seq = (tmp_seq << 4) + 4;
                            }
                        }
                        // At right bound
                        if (((contig_offset + atol(snp_position) - 1) / 8 + (SEQ_CHUNK_NUMBER / 2) + i + 1) * 8 < contig_offset + mi->seq[seq_num].len &&
                            ((contig_offset + atol(snp_position) - 1) / 8 + (SEQ_CHUNK_NUMBER / 2) + i + 2) * 8 > contig_offset + mi->seq[seq_num].len) {
                            tmp_seq = tmp_seq << (4 * (8 - (contig_offset + mi->seq[seq_num].len) % 8));
                            for (int j = 0; j < 8 - (contig_offset + mi->seq[seq_num].len) % 8; j++)
                                tmp_seq = (tmp_seq >> 4) | 1073741824; // FIRST N
                        }
                    }
                }

                for (int j = 0; j < 8; j++) {
                    uint32_t tmp = tmp_seq % 16;
                    switch (tmp) {
                        case 0:
                            original_ref_seq_ext[i * 8 + j] = 'A';
                            break;
                        case 1:
                            original_ref_seq_ext[i * 8 + j] = 'C';
                            break;
                        case 2:
                            original_ref_seq_ext[i * 8 + j] = 'G';
                            break;
                        case 3:
                            original_ref_seq_ext[i * 8 + j] = 'T';
                            break;
                        case 4:
                            original_ref_seq_ext[i * 8 + j] = 'N';
                    }
                    tmp_seq = tmp_seq / 16;
                }
            }

            //Create new window
            char * new_ref_seq;
            new_ref_seq = malloc(sizeof(char) * (SEQ_CHUNK_NUMBER * 8 + 1 - strlen(str[3]) + 1 + ext_chunk_count * 8));
            memcpy(new_ref_seq, original_ref_seq, EXTRA_GAP + SIDE_SIZE + (contig_offset + atol(snp_position) - 1) % 8);
            new_ref_seq[EXTRA_GAP + SIDE_SIZE + (contig_offset + atol(snp_position) - 1) % 8] = snp_to;
            new_ref_seq[EXTRA_GAP + SIDE_SIZE + (contig_offset + atol(snp_position) - 1) % 8 + 1] = '\0';
            new_ref_seq = strcat(new_ref_seq, &original_ref_seq[EXTRA_GAP + SIDE_SIZE + (contig_offset + atol(snp_position) - 1) % 8 + strlen(str[3])]);
            new_ref_seq = strcat(new_ref_seq, original_ref_seq_ext);

            //Finds minimizer in window
            mm128_v minimizer_array = {0, 0, 0};
            mm_sketch(0, &new_ref_seq[(contig_offset + atol(snp_position) - 1) % 8], SIDE_SIZE * 2 + 1, mi->w, mi->k,
                      0, mi->flag & MM_I_HPC, &minimizer_array);


            free(original_ref_seq_ext);
            free(new_ref_seq);

            for (int i = 0; i < minimizer_array.n; i++) {
                if (minimizer_array.a[i].y < SIDE_SIZE * 2 + 2) continue; // TODO change 50 to 48 for similarity with SNPS (it will take extra calculation but no changes)
                if (minimizer_array.a[i].y > (SIDE_SIZE + mi->k) * 2 - 1) continue;
                minimizer_array.a[i].y = (seq_num << 32) + (atol(snp_position) - SIDE_SIZE - 1 + minimizer_array.a[i].y / 2) * 2 +
                                         (minimizer_array.a[i].y % 2) + (strlen(str[3]) - 1) * 2;
                mm_idx_push(mi, minimizer_array.a[i].x, minimizer_array.a[i].y);
            }
        } else if (strlen(str[3]) == 1 && strlen(str[4]) > 1) {
            //Insertion
            //Create new window
            char * new_ref_seq;
            new_ref_seq = malloc(sizeof(char) * (SEQ_CHUNK_NUMBER * 8 + 1 + strlen(str[4]) - 1));
            memcpy(new_ref_seq, original_ref_seq, EXTRA_GAP + SIDE_SIZE + (contig_offset + atol(snp_position) - 1) % 8);
            new_ref_seq[EXTRA_GAP + SIDE_SIZE + (contig_offset + atol(snp_position) - 1) % 8] = '\0';
            new_ref_seq = strcat(new_ref_seq, str[4]);
            new_ref_seq = strcat(new_ref_seq, &original_ref_seq[EXTRA_GAP + SIDE_SIZE + (contig_offset + atol(snp_position) - 1) % 8 + 1]);

            //Finds minimizer in window
            mm128_v minimizer_array = {0, 0, 0};
            mm_sketch(0, &new_ref_seq[(contig_offset + atol(snp_position) - 1) % 8], SIDE_SIZE * 2 + 1 + strlen(str[4]) - 1, mi->w, mi->k,
                      0, mi->flag & MM_I_HPC, &minimizer_array);

            for (int i = 0; i < minimizer_array.n; i++) {
                if (minimizer_array.a[i].y < SIDE_SIZE * 2 + 2) continue;
                if (minimizer_array.a[i].y > (SIDE_SIZE + mi->k) * 2 - 1 + (strlen(str[4]) - 1) * 2) continue;
                minimizer_array.a[i].y = (seq_num << 32) + (atol(snp_position) - SIDE_SIZE - 1 + minimizer_array.a[i].y / 2) * 2 +
                                         (minimizer_array.a[i].y % 2);
                mm_idx_push(mi, minimizer_array.a[i].x, minimizer_array.a[i].y);
            }
            free(new_ref_seq);
        } else {
            printf("This mutation doesn't support yet.\n");
        }
    }
    fclose(vcf);
}
