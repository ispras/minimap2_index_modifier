#ifndef MM_VARALT_H
#define MM_VARALT_H

/* --var-alt: re-score a candidate's mismatches against known panel SNV ALTs.
 * A mismatch whose read base is a known ALT at that base position counts as a
 * match in dp_score/dp_max/mlen; CIGAR, positions and NM are unchanged.
 *
 * Table (plain or gzip): <seq name> TAB <1-based pos> TAB <ALTs over ACGT>.
 * ALTs equal to the base are dropped; duplicate positions are merged. */

#include <stdint.h>
#include "minimap.h"

typedef struct {
	uint32_t n;       /* positions on this sequence */
	uint32_t *pos;    /* 0-based, strictly increasing */
	uint8_t *alt;     /* bit c (c = 0..3 for A,C,G,T) set: c is a known ALT at pos[i] */
} mm_va_seq_t;

typedef struct mm_varalt_s {
	uint32_t n_seq;           /* == mi->n_seq */
	mm_va_seq_t *s;
	uint32_t *pos_pool;       /* one allocation for all sequences */
	uint8_t *alt_pool;
	int64_t n_pos;            /* positions kept */
	int64_t n_alleles;        /* ALT alleles kept (bits set) */
	int64_t n_lines;          /* data lines read */
	int64_t n_unknown_seq;    /* lines whose sequence is not in the index */
	int64_t n_bad;            /* malformed lines (skipped) */
	int64_t n_eq_ref;         /* ALT letters equal to the base at that position (dropped) */
	int64_t n_out_of_range;   /* positions beyond the sequence end (skipped) */
	int64_t n_dup;            /* lines merged into an earlier line of the same position */
} mm_varalt_t;

#ifdef __cplusplus
extern "C" {
#endif

/* Load FILE against the index; NULL on I/O error. */
mm_varalt_t *mm_varalt_load(const mm_idx_t *mi, const char *fn);
void mm_varalt_destroy(mm_varalt_t *va);

/* Known-ALT mask at (rid, 0-based pos); 0 if none. */
int mm_varalt_mask(const mm_varalt_t *va, int32_t rid, uint32_t pos);

/* Re-score one candidate. qstr: the forward read; the CIGAR must be against
 * base coordinates. Sets r->va and returns it (0: nothing changed). */
int mm_varalt_rescore(const mm_idx_t *mi, const mm_mapopt_t *opt, int qlen, const char *qstr, mm_reg1_t *r);

#ifdef __cplusplus
}
#endif

#endif
