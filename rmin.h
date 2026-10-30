/* rmin: filter the minimizers added by the index modifier (add_variants) by the
 * minimizers/k-mers seen in the reads. Base reference minimizers are never
 * touched; RMIN_NONE leaves the index unchanged. Count file format: rmin_count.c. */
#ifndef RMIN_H
#define RMIN_H
#include <stdint.h>
#include <stdio.h>
#include "minimap.h"

#define RMIN_NONE   0
#define RMIN_STRICT 1  /* keep a minimizer whose hash has read-minimizer count >= C */
#define RMIN_COMBO  2  /* keep a whole combination if >= F of its minimizers have count >= C */
#define RMIN_KMER   3  /* as strict, against all canonical read k-mers */
#define RMIN_RAND   4  /* drop a minimizer with probability R (hash of x, y, seed) */
#define RMIN_LIST   5  /* keep a combination if all its SNVs are in a CHR POS ALT list */

typedef struct {
	int mode, min_cnt;
	double frac, rate;
	uint64_t seed;
	const char *fn_mm, *fn_km, *fn_list, *fn_dump;
	const char *fn_ref; /* --read-min-ref-kmers: reference k-mers (w=1); an ALT k-mer found
	                       there is dropped in any mode */
} rmin_opt_t;

extern rmin_opt_t rmin_opt;

int  rmin_parse_mode(const char *s);
int  rmin_init(void);           /* load files; 0 on success */
int  rmin_active(void);
void rmin_rec_begin(const char *chr, unsigned long curr_pos);
void rmin_rec_end(void);
/* a[0..n): candidates of one add_variants() combination (y in sequence
 * coordinates); sets keep[i] to 1/0 and returns the number kept. */
int  rmin_filter(const mm128_t *a, int n, uint8_t *keep, const char *chr, unsigned long curr_pos,
                 char **REF_arr, char **ALT_arr, unsigned long *POS_all, int N_SNP, int k);
void rmin_report(FILE *fp);
void rmin_destroy(void);

#endif
