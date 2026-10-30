#ifndef MM_HAPDP_H
#define MM_HAPDP_H

/* --hap-dp: haplotype-consistent DP score. Within a candidate's window, each
 * distinct panel-haplotype allele string over known SNVs is substituted into
 * the target and scored; the best string sets dp_max/dp_score/mlen. CIGAR,
 * positions and NM stay against the base reference.
 *
 * Table (plain or gzip), one line per (position, ALT):
 *     <seq name> TAB <1-based pos> TAB <ALT: ACGT> TAB <mask_hi> TAB <mask_lo>
 * masks are hex carrier sets (bit h of lo = haplotype h < 64, bit h-64 of hi);
 * '#' lines are comments, "#n_hap<TAB>N" declares the haplotype count.
 * --hap-dp-mask FILE lists 0-based haplotype indices allowed to contribute. */

#include <stdint.h>
#include "minimap.h"

typedef struct {
	uint32_t pos;     /* 0-based */
	uint32_t alt;     /* 0..3 = A,C,G,T */
	uint64_t m[2];    /* carriers: bit h of m[h>>6] */
} mm_hd_rec_t;

typedef struct {
	uint32_t n;
	mm_hd_rec_t *a;   /* sorted by pos, then alt */
} mm_hd_seq_t;

typedef struct mm_hapdp_s {
	uint32_t n_seq;
	mm_hd_seq_t *s;
	mm_hd_rec_t *pool;
	int32_t n_hap;            /* declared (#n_hap) or inferred */
	int32_t max_str;          /* max strings per window including the base */
	int32_t has_allow;
	uint64_t allow[2];        /* --hap-dp-mask */
	int64_t n_rec;            /* records kept */
	int64_t n_lines, n_bad, n_unknown_seq, n_out_of_range, n_eq_ref, n_nocarrier, n_dup, n_conflict;
} mm_hapdp_t;

/* one window: n_rec records, n_str strings (the base string is NOT among them);
 * on[s * n_rec + j] = 1 iff string s carries the ALT of record j */
typedef struct {
	int32_t n_rec, n_str, n_hap_win;  /* n_hap_win: haplotypes carrying >= 1 ALT of the window */
	int32_t truncated;                /* strings dropped by max_str */
	const mm_hd_rec_t *rec;
	uint8_t *on;
	int32_t *cnt;                     /* haplotypes per string */
} mm_hd_win_t;

#ifdef __cplusplus
extern "C" {
#endif

mm_hapdp_t *mm_hapdp_load(const mm_idx_t *mi, const char *fn);
int mm_hapdp_load_mask(mm_hapdp_t *hd, const char *fn);   /* 0 on success */
void mm_hapdp_destroy(mm_hapdp_t *hd);

/* records of rid in [st, en): index of the first and count; 0 = none (the common case) */
int32_t mm_hapdp_range(const mm_hapdp_t *hd, int32_t rid, int32_t st, int32_t en, uint32_t *first);

/* enumerate the distinct haplotype strings over the records of [st, en); buffers from km
 * (free with mm_hapdp_win_free).  Returns n_str (0: no string differs from the base). */
int32_t mm_hapdp_window(void *km, const mm_hapdp_t *hd, int32_t rid, int32_t st, int32_t en, mm_hd_win_t *w);
void mm_hapdp_win_free(void *km, mm_hd_win_t *w);

/* best ungapped score over the base and the strings, for qseq/tseq of [st, st+len) on rid
 * (tseq = base).  Used for the ungapped-vs-ksw decision of the short-read middle. */
int32_t mm_hapdp_ungapped_best(void *km, const mm_hapdp_t *hd, const mm_mapopt_t *opt, int32_t rid, int32_t st, int32_t len,
							   const uint8_t *qseq, const uint8_t *tseq, int32_t base_score);

/* Re-score a finished candidate (CIGAR against the base, after mm_update_extra()).
 * qseq: the read in alignment orientation starting at the aligned query start;
 * mat: the 5x5 scoring matrix; log_gap as in mm_update_extra().  Sets r->hd,
 * r->hd_dmlen and updates dp_max / dp_score / mlen if a string beats the base.
 * Returns 1 if the candidate changed. */
int mm_hapdp_rescore(void *km, const mm_idx_t *mi, const mm_mapopt_t *opt, const int8_t *mat, int log_gap,
					 const uint8_t *qseq, mm_reg1_t *r);

/* statistics: an ungapped middle accepted only thanks to a haplotype string */
void mm_hapdp_count_mid_flip(void);

#ifdef __cplusplus
}
#endif

#endif
