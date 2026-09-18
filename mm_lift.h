#ifndef MM_LIFT_H
#define MM_LIFT_H

/* ---------------------------------------------------------------------------
 * In-aligner projection of augmentation-contig coordinates onto the base
 * reference.
 *
 * WHY THIS EXISTS
 * ---------------
 * The augmented reference is T' = T + C, where T is GRCh38 (195 sequences) and
 * C is 569,870 short contigs, each carrying one known indel's ALT allele
 * between 100 bp of its native flanks.  Today a read that lands on a member of
 * C is rewritten into T coordinates *after* alignment by levioSAM2, driven by
 * a separate chain file.  That costs 11 h 09 m per 40x genome
 * (scratchpad/remote/costs.md, table b.1, sourced from bgi/log/progress.txt)
 * and the chain is an independent artifact that can — and does — disagree with
 * the FASTA: only 195 of 570,065 chain records are internally consistent
 * (0.0342 %, log/t1pilot/chainfix.log step [3]), and those 195 are the
 * whole-chromosome pass-through records, i.e. *every* contig record is wrong.
 *
 * This header defines the projection as data owned by the index, so that
 * minimap2 can emit base-reference coordinates directly and no second artifact
 * exists to drift.
 *
 * THE GEOMETRY, DERIVED FROM THE GENERATOR
 * ----------------------------------------
 * gen2/contig_generation.FIXED.py builds each contig from a VCF record
 * (chrom, pos, REF, ALT) with contig_half = 100 as:
 *
 *     contig = T[pos-100 .. pos-1] ++ ALT ++ T[pos+|REF| .. pos+|REF|+99]
 *              \__ 100 bp left __/         \______ 100 bp right flank ____/
 *
 * VCF indel convention puts the anchor base at pos in both REF[0] and ALT[0],
 * so contig[0..100] is 101 contiguous bases of T ending at pos.  Therefore the
 * projection is a two-block affine map with a single junction:
 *
 *     block A : contig [0, jct)                -> T0 [off_a,            off_a + jct)
 *     (gap)   : contig [jct, jct+ins_len)      -> nothing   (novel ALT bases)
 *               T0     [off_a+jct, +del_len)   -> nothing   (deleted REF bases)
 *     block B : contig [jct+ins_len, len)      -> T0 [off_a + jct + del_len, ...)
 *
 * with, for this reference build,
 *     jct     = 101                 (contig_half + 1)
 *     off_a   = pos - 101           (0-based T coordinate of contig offset 0)
 *     ins_len = max(|ALT| - |REF|, 0)
 *     del_len = max(|REF| - |ALT|, 0)
 *
 * Exactly one of ins_len / del_len is non-zero: the generator's filter keeps
 * only records where |ALT| > |REF| (insertion) or |ALT| < |REF| (deletion).
 *
 * Worked against the real chain (data/chain_sample.txt, from gen2/fixed2.chain):
 *   chr1_contig_30884, len 205, chain "101 4<TAB>0" -> insertion, ins_len 4,
 *      off_a = 30884-101 = 30783, which is the chain's own qStart.
 *   chr1_contig_54750, len 201, chain "101 0<TAB>4" -> deletion, del_len 4,
 *      off_a = 54750-101 = 54649, again the chain's qStart.
 * The chain's *block lengths* are wrong (99 instead of 100 for insertions,
 * 98+|REF| instead of 100 for deletions) but its header coordinates are right,
 * which is why the map below is derivable and the chain is not trustworthy.
 *
 * WHAT IS AND IS NOT DERIVABLE FROM THE INDEX AS IT STANDS
 * --------------------------------------------------------
 * mm_idx_dump() already stores every sequence's name and length
 * (index.c:598-607), so (name, len) is free at load time.  From those alone:
 *     pos     = atoi(strstr(name, "_contig_") + 8)      -> off_a = pos - 101
 *     |ALT|   = len - 200
 * and for an insertion |REF| = 1, so ins_len = len - 201.  That closes the map
 * for every insertion contig.  It does NOT close it for a deletion, because a
 * deletion's |REF| is nowhere in (name, len).
 *
 * Measured on the production reference (data/fixed2.fasta.fai, 569,870
 * contigs), the ambiguity is almost empty:
 *     len == 201                348,334 contigs  -> |ALT| = 1, must be a
 *                                                   deletion (the generator
 *                                                   requires |ALT| >= 5 for
 *                                                   insertions), |REF| unknown
 *     len in {203, 204}                  91      -> same, |REF| unknown
 *     len >= 205                   221,445      -> 221,283 insertions (fully
 *                                                   determined) + 162 deletions
 * Total deletions 348,587 and insertions 221,283, which reproduces
 * log/t1pilot/result.txt exactly.  So del_len must be supplied for 348,587
 * contigs and only 162 of them are length-ambiguous with an insertion.
 *
 * CONCLUSION: one int32 per contig (the signed indel delta d = |ALT| - |REF|)
 * is the whole missing payload.  At 569,870 contigs that is 2.28 MB, against a
 * 9,479,433,012-byte index (costs.md table c) -- 0.024 %.  Compare the .clft
 * chain index it replaces: 271,714,958 bytes, built in 82 s
 * (log/t1pilot/chainfix.log step [4]).
 * ------------------------------------------------------------------------- */

#include <stdint.h>
#include "minimap.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Per-sequence projection record.  One entry per sequence in mm_idx_t::seq,
 * indexed by rid, so lookup is an array subscript on the hot path. */
typedef struct {
	int32_t base_rid;  /* rid of the base sequence this contig projects onto;
	                    * -1 means "this sequence IS a base sequence" and the
	                    * projection is the identity. */
	int32_t jct;       /* contig offset of the first base after block A.
	                    * 101 for this reference build.  Meaningless when
	                    * base_rid < 0. */
	int32_t ins_len;   /* contig bases at [jct, jct+ins_len) with no counterpart
	                    * in the base reference (novel ALT bases). >= 0 */
	int32_t del_len;   /* base-reference bases at [off_a+jct, +del_len) with no
	                    * counterpart in the contig (deleted REF bases). >= 0 */
	int64_t off_a;     /* 0-based base-reference coordinate of contig offset 0 */
} mm_lift_rec_t;

/* Projection table.  Owned by the index; n == mi->n_seq, indexed by rid. */
typedef struct mm_lift_s {
	int32_t n;
	mm_lift_rec_t *a;
} mm_lift_t;

/* off_b: base-reference coordinate of contig offset 0 *for block B*.
 * Derived, never stored, so the two cannot drift.
 *    contig c in block B  ->  off_b(rec) + c
 * Identity: off_b = off_a + del_len - ins_len. */
static inline int64_t mm_lift_off_b(const mm_lift_rec_t *rec)
{
	return rec->off_a + rec->del_len - rec->ins_len;
}

/* Reference-consuming span of a packed CIGAR (M/D/N/=/X).  Used to recompute
 * mm_reg1_t::re after the splice, so that re always equals rs + span and the
 * two cannot drift. */
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

/* Zone of a contig offset.  0 = block A, 1 = the novel-ALT gap, 2 = block B. */
#define MM_LIFT_ZONE_A   0
#define MM_LIFT_ZONE_INS 1
#define MM_LIFT_ZONE_B   2

static inline int mm_lift_zone(const mm_lift_rec_t *rec, int64_t c)
{
	if (c < rec->jct) return MM_LIFT_ZONE_A;
	if (c < (int64_t)rec->jct + rec->ins_len) return MM_LIFT_ZONE_INS;
	return MM_LIFT_ZONE_B;
}

/* ---------------------------------------------------------------------------
 * Name parsing.
 *
 * Contract: returns 1 and fills *chrom_len (length of the chromosome-name
 * prefix within name) and *pos (the 1-based variant coordinate) iff name
 * matches "<chrom>_contig_<digits>" with at least one digit and nothing after
 * them.  Returns 0 otherwise, leaving the outputs untouched.
 *
 * Uses the LAST occurrence of "_contig_", so a chromosome whose own name
 * contains the token cannot shadow the split.  strstr()-based parsing in
 * map.c:256 uses the first occurrence; that is a latent difference and this
 * function is the one to keep.
 *
 * NOTE this deliberately does NOT reproduce the truncation bug in
 * gen2/contig_generation.py (the superseded generator), whose
 * `contig_header[1:5]` turns ">chr10_contig_N" into "chr1".  The production
 * generator contig_generation.FIXED.py splits on "_" instead and is correct;
 * the bug is recorded here only so a reader comparing the two is not confused.
 * ------------------------------------------------------------------------- */
int mm_lift_parse_name(const char *name, int32_t *chrom_len, int64_t *pos);

/* ---------------------------------------------------------------------------
 * Project a single contig offset onto the base reference.
 *
 * Contract:
 *   rec must have base_rid >= 0.
 *   c is a 0-based contig offset in [0, contig_len].
 *   Returns the 0-based base-reference coordinate, choosing the block that
 *   contains c.  For c in the novel-ALT gap there is no counterpart; the
 *   convention is to return the junction's base-reference coordinate
 *   (off_a + jct), i.e. to collapse the gap to a point.  Callers that must
 *   distinguish should test mm_lift_zone() first.
 * ------------------------------------------------------------------------- */
int64_t mm_lift_pos(const mm_lift_rec_t *rec, int64_t c);

/* ---------------------------------------------------------------------------
 * Project an alignment: rewrite start coordinate and CIGAR from contig space
 * into base-reference space.
 *
 * Input:
 *   rec        projection record, base_rid >= 0
 *   ctg_len    length of the contig (mi->seq[rid].len)
 *   rs         0-based contig start of the alignment (mm_reg1_t::rs)
 *   cigar,n_cigar   packed minimap2/BAM CIGAR (len<<4 | op), as in
 *                   mm_extra_t::cigar.  Soft/hard clips are NOT included:
 *                   minimap2 keeps clipping outside r->p->cigar and synthesises
 *                   it in write_sam_cigar() (format.c:489-491) from r->qs/r->qe.
 *   qs, qe     query interval of the alignment (mm_reg1_t::qs/qe)
 *
 * Output (all written only on success):
 *   *rs_out        0-based base-reference start
 *   *qs_out,*qe_out  possibly tightened query interval: when the alignment
 *                  begins or ends inside the novel-ALT gap those query bases
 *                  have nowhere to go and are converted to clipping by moving
 *                  qs/qe.  Callers must write these back into the mm_reg1_t so
 *                  that write_sam_cigar()'s clip_len[] stays consistent.
 *   cigar_out      caller-allocated buffer, capacity cap ops
 *   *n_out         number of ops written
 *
 * Returns:
 *    0  success
 *   -1  buffer too small (cap insufficient); *n_out holds the required count
 *   -2  the alignment lies entirely inside the novel-ALT gap, so it has no
 *       base-reference representation at all.  The caller must decide policy
 *       (drop, emit unmapped, or retain contig coordinates).  Cannot occur for
 *       this reference with >=150 bp reads: the generator rejects ALT longer
 *       than 100 bp (contig_generation.FIXED.py, `if len(fields[4]) > 100`), so
 *       the gap is at most 99 bp and cannot contain a whole read.
 *
 * Guarantees on the emitted CIGAR:
 *   G1  reference-consuming length equals the base-reference span implied by
 *       the projection (no invented and no lost reference bases);
 *   G2  query-consuming length equals qe_out - qs_out;
 *   G3  no zero-length ops, no two adjacent ops with the same opcode;
 *   G4  does not begin or end with D or N;
 *   G5  a `del_len D` op appears at the junction if and only if the alignment
 *       has reference-consuming coverage strictly on both sides of it.  This is
 *       the op levioSAM2 emits unconditionally from the chain, which is the
 *       measured +L pathology when the FASTA and the chain disagree.
 *
 * NOT guaranteed, and deliberately so: indel placement is not normalised
 * against the base reference.  See DESIGN.md class C3 -- left-shifting a
 * junction indel to its leftmost equivalent position is a separate, explicit
 * pass that needs mm_idx_t::S, and it is not folded in here.
 * ------------------------------------------------------------------------- */
int mm_lift_cigar(const mm_lift_rec_t *rec, int32_t ctg_len,
                  int64_t rs, const uint32_t *cigar, int32_t n_cigar,
                  int32_t qs, int32_t qe,
                  int64_t *rs_out, int32_t *qs_out, int32_t *qe_out,
                  uint32_t *cigar_out, int32_t cap, int32_t *n_out);

/* ---------------------------------------------------------------------------
 * Table construction.
 *
 * mm_lift_build(mi, fn_map, jct) fills mi->lift with one record per sequence.
 *
 *   - A sequence whose name does not parse as "<chrom>_contig_<pos>" gets
 *     base_rid = -1 (identity).
 *   - A contig whose chromosome prefix does not name a sequence in the index
 *     also gets base_rid = -1, and is counted as unresolved: it cannot be
 *     projected, so it is left alone rather than projected wrongly.
 *   - off_a = pos - jct, jct as given (101 for this reference build).
 *   - ins_len/del_len come from the signed delta d = |ALT| - |REF|:
 *       d > 0 -> ins_len = d,  del_len = 0
 *       d < 0 -> ins_len = 0,  del_len = -d
 *     d is read from fn_map, a two-column "<contig name>\t<d>" TSV.  For a
 *     contig absent from the file, d is inferred as len-201 when that is >= 4
 *     (an insertion, since the generator requires |ALT| >= 5 and |REF| == 1);
 *     otherwise the contig is unresolved.  See the metadata audit above for
 *     why inference alone is not enough: every deletion contig needs the file.
 *
 * Returns the number of sequences with base_rid >= 0, or -1 on error.
 * Writes the unresolved count to *n_unresolved when that is non-NULL.
 * ------------------------------------------------------------------------- */
int mm_lift_build(mm_idx_t *mi, const char *fn_map, int32_t jct, int64_t *n_unresolved);

void mm_lift_destroy(mm_lift_t *lt);

/* ---------------------------------------------------------------------------
 * Project a whole candidate set in place.
 *
 * For every reg whose rid carries a projection record, rewrites rid, rs, re,
 * qs, qe and the CIGAR in r->p, and sets r->lifted.  Regs on base sequences
 * and regs whose projection fails are left exactly as they were, so the
 * function is safe to call unconditionally.
 *
 * Must be called AFTER mapping quality is set -- mapq depends on the candidate
 * set as the aligner scored it, and after projection the name-based rule in
 * collect_seed_chromosome_names() can no longer tell a contig from a
 * chromosome (every reg is on a chromosome by then).  Read r->lifted instead.
 *
 * Must be called BEFORE mm_pair(), so that a read on a contig and its mate on
 * the chromosome can be judged a proper pair: before projection they carry
 * different rids and mm_pair() cannot relate them at all.
 *
 * Returns the number of regs actually projected.
 * ------------------------------------------------------------------------- */
int mm_lift_regs(const mm_idx_t *mi, int n_regs, mm_reg1_t *regs);

#ifdef __cplusplus
}
#endif

#endif /* MM_LIFT_H */
