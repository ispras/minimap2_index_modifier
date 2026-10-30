#ifndef MM_LIFT_H
#define MM_LIFT_H

/* Projection of augmentation-contig coordinates onto the base reference.
 *
 * A contig is a base-reference window with one or more indel ALT alleles
 * spliced in. Its map to the base is piecewise affine: blocks separated by
 * gaps, each gap inserting ins_len novel bases and/or skipping del_len bases. */

#include <stdint.h>
#include "minimap.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Gap g is (jct_g, ins_g, del_g) in contig coordinates; jct_g is the contig
 * offset of the first base after block g. Block offsets are derived:
 * off_0 = off_a, off_{b+1} = off_b + del_b - ins_b. Gap 0 is stored in the
 * record, further gaps in the pool mm_lift_t::xg. */
#define MM_LIFT_MAX_GAPS 8   /* gaps per contig, first included; the build rejects more */

typedef struct {
	int32_t jct, ins_len, del_len;   /* one junction, contig coordinates, as in mm_lift_rec_t */
} mm_lift_gap_t;

/* Per-sequence projection record, indexed by rid */
typedef struct {
	int32_t base_rid;  /* base sequence rid; -1 = this is a base sequence (identity) */
	int32_t jct;       /* contig offset of the first base after block A */
	int32_t ins_len;   /* novel contig bases at [jct, jct+ins_len) */
	int32_t del_len;   /* base bases at [off_a+jct, +del_len) absent from the contig */
	int64_t off_a;     /* 0-based base-reference coordinate of contig offset 0 */
	int32_t n_xgap;    /* number of gaps after the first */
	const mm_lift_gap_t *xgap;   /* ascending jct; points into mm_lift_t::xg, NULL if n_xgap == 0 */
} mm_lift_rec_t;

/* Projection table owned by the index; n == mi->n_seq */
typedef struct mm_lift_s {
	int32_t n;
	mm_lift_rec_t *a;
	int64_t n_jct_from_map;   /* contigs whose jct came from the map file */
	mm_lift_gap_t *xg;        /* pool of the extra gaps */
	int64_t n_xg;             /* entries in xg */
	int64_t n_multi;          /* contigs with more than one gap */
} mm_lift_t;

/* base coordinate of contig offset 0 for block B (after the first gap) */
static inline int64_t mm_lift_off_b(const mm_lift_rec_t *rec)
{
	return rec->off_a + rec->del_len - rec->ins_len;
}

/* number of gaps; 1 for a single-indel contig */
static inline int32_t mm_lift_n_gaps(const mm_lift_rec_t *rec)
{
	return 1 + (rec->xgap != NULL ? rec->n_xgap : 0);
}

/* gap g (0-based) */
static inline mm_lift_gap_t mm_lift_gap_at(const mm_lift_rec_t *rec, int32_t g)
{
	mm_lift_gap_t x;
	if (g == 0) { x.jct = rec->jct; x.ins_len = rec->ins_len; x.del_len = rec->del_len; }
	else x = rec->xgap[g - 1];
	return x;
}

/* base coordinate of contig offset 0 for block b (0..n_gaps) */
static inline int64_t mm_lift_block_off(const mm_lift_rec_t *rec, int32_t b)
{
	int64_t o = rec->off_a; int32_t g;
	for (g = 0; g < b; ++g) { mm_lift_gap_t x = mm_lift_gap_at(rec, g); o += x.del_len - x.ins_len; }
	return o;
}

/* reference-consuming span of a packed CIGAR (M/D/N/=/X) */
static inline int64_t mm_lift_ref_span(const uint32_t *cigar, int32_t n)
{
	int64_t s = 0; int32_t i;
	for (i = 0; i < n; ++i) {
		int op = cigar[i] & 0xf;
		if (op == 0 || op == 2 || op == 3 || op == 7 || op == 8)
			s += cigar[i] >> 4;
	}
	return s;
}

/* query-consuming span of a packed CIGAR (M/I/=/X) */
static inline int64_t mm_lift_query_span(const uint32_t *cigar, int32_t n)
{
	int64_t s = 0; int32_t i;
	for (i = 0; i < n; ++i) {
		int op = cigar[i] & 0xf;
		if (op == 0 || op == 1 || op == 7 || op == 8)
			s += cigar[i] >> 4;
	}
	return s;
}

/* Zone of a contig offset: 2b = block b, 2g+1 = novel bases of gap g */
#define MM_LIFT_ZONE_A   0
#define MM_LIFT_ZONE_INS 1
#define MM_LIFT_ZONE_B   2

static inline int mm_lift_zone(const mm_lift_rec_t *rec, int64_t c)
{
	int32_t g, n = mm_lift_n_gaps(rec);
	for (g = 0; g < n; ++g) {
		mm_lift_gap_t x = mm_lift_gap_at(rec, g);
		if (c < x.jct) return 2 * g;
		if (c < (int64_t)x.jct + x.ins_len) return 2 * g + 1;
	}
	return 2 * n;
}

/* first contig offset AFTER zone z; INT64_MAX for the last block */
static inline int64_t mm_lift_zone_end(const mm_lift_rec_t *rec, int z)
{
	int32_t n = mm_lift_n_gaps(rec);
	mm_lift_gap_t x;
	if (z >= 2 * n) return INT64_MAX;
	x = mm_lift_gap_at(rec, z / 2);
	return (z & 1) ? (int64_t)x.jct + x.ins_len : (int64_t)x.jct;
}

/* Parse "<chrom>_contig_<pos>" or "<chrom>_contig_<pos>_h<serial>" (last
 * "_contig_" wins). On match returns 1 and sets *chrom_len (prefix length) and
 * *pos (1-based); otherwise returns 0 and leaves the outputs untouched. */
int mm_lift_parse_name(const char *name, int32_t *chrom_len, int64_t *pos);

/* Project 0-based contig offset c (rec->base_rid >= 0) to a 0-based base
 * coordinate. An offset inside a novel-ALT gap maps to the gap's junction. */
int64_t mm_lift_pos(const mm_lift_rec_t *rec, int64_t c);

/* Project an alignment: rewrite start and CIGAR from contig to base space.
 *
 * rs is the 0-based contig start; cigar is packed (len<<4|op) without clips;
 * qs/qe/rev are as in mm_reg1_t (the CIGAR is in reference order, qs/qe in
 * forward-read order). Query bases at either end that fall in a novel-ALT gap
 * become clipping via *qs_out and *qe_out, which the caller must write back.
 * cigar_out is caller-allocated with capacity cap.
 *
 * Returns 0 on success, -1 if cap is too small (*n_out = required count), -2 if
 * the alignment lies entirely inside a novel-ALT gap. Outputs are written only
 * on success. The emitted CIGAR has no zero-length or repeated adjacent ops,
 * does not start or end with D/N, and contains a gap's deletion only when it is
 * covered on both sides. Indels are not left-normalised here. */
int mm_lift_cigar(const mm_lift_rec_t *rec, int32_t ctg_len,
                  int64_t rs, const uint32_t *cigar, int32_t n_cigar,
                  int32_t qs, int32_t qe, int rev,
                  int64_t *rs_out, int32_t *qs_out, int32_t *qe_out,
                  uint32_t *cigar_out, int32_t cap, int32_t *n_out);

/* Build mi->lift from contig names and the map file fn_map (TSV):
 *   <name> <d>  or  <name> <ins_len> <del_len> [<jct> [<jct:ins:del;...>]]
 * d = |ALT|-|REF|; the optional 4th column is the junction (default jct), the
 * 5th lists further gaps ascending (at most MM_LIFT_MAX_GAPS in total). A
 * contig absent from the file is an insertion of len-201 if that is >= 4.
 * Sequences that do not parse, or name an unknown chromosome, or have invalid
 * geometry, get base_rid = -1 (counted in *n_unresolved if non-NULL).
 * Returns the number of projectable sequences, or -1 on error. */
int mm_lift_build(mm_idx_t *mi, const char *fn_map, int32_t jct, int64_t *n_unresolved);

void mm_lift_destroy(mm_lift_t *lt);

/* Project all regs on contigs in place (rid, rs, re, qs, qe, CIGAR), setting
 * r->lifted and r->lift_rid/r->lift_rs. Regs on base sequences, or whose
 * projection fails, are left unchanged. Call before the mapq rule and
 * mm_pair(). Returns the number of regs projected. */
int mm_lift_regs(const mm_idx_t *mi, int n_regs, mm_reg1_t *regs);

/* Return 1 if a and b (after projection) are on the same sequence and strand
 * and share at least one (query, reference) base pair, i.e. a common diagonal;
 * uses the CIGARs when present, otherwise the end points. */
int mm_lift_same_locus(const mm_reg1_t *a, const mm_reg1_t *b);

#ifdef __cplusplus
}
#endif

#endif /* MM_LIFT_H */
