/* Projection of augmentation-contig coordinates onto the base reference */

#include <string.h>
#include <stdlib.h>
#include "mm_lift.h"

/* CIGAR opcodes, as MM_CIGAR_* in minimap.h */
#define OP_M 0
#define OP_I 1
#define OP_D 2
#define OP_N 3
#define OP_S 4
#define OP_H 5
#define OP_P 6
#define OP_EQ 7
#define OP_X 8

#define CIG_LEN(c) ((c) >> 4)
#define CIG_OP(c)  ((c) & 0xf)
#define CIG_MK(l, o) ((uint32_t)(l) << 4 | (o))

static inline int op_consumes_ref(int op)
{
	return op == OP_M || op == OP_D || op == OP_N || op == OP_EQ || op == OP_X;
}

static inline int op_consumes_query(int op)
{
	return op == OP_M || op == OP_I || op == OP_S || op == OP_EQ || op == OP_X;
}

int mm_lift_parse_name(const char *name, int32_t *chrom_len, int64_t *pos)
{
	const char *tag = NULL, *p, *q;
	static const char TOK[] = "_contig_";
	const size_t TOKLEN = sizeof(TOK) - 1;
	int64_t v;

	if (name == NULL) return 0;
	/* last occurrence */
	for (p = name; (q = strstr(p, TOK)) != NULL; p = q + 1) tag = q;
	if (tag == NULL || tag == name) return 0;  /* no token, or empty chrom */

	p = tag + TOKLEN;
	if (*p < '0' || *p > '9') return 0;        /* need at least one digit */
	v = 0;
	for (; *p >= '0' && *p <= '9'; ++p) {
		if (v > (INT64_MAX - (*p - '0')) / 10) return 0;  /* overflow */
		v = v * 10 + (*p - '0');
	}
	if (p[0] == '_' && p[1] == 'h') {          /* "_h<digits>" suffix */
		const char *s = p + 2;
		if (*s < '0' || *s > '9') return 0;
		while (*s >= '0' && *s <= '9') ++s;
		if (*s != '\0') return 0;
		p = s;
	}
	if (*p != '\0') return 0;                  /* trailing junk */

	*chrom_len = (int32_t)(tag - name);
	*pos = v;
	return 1;
}

int64_t mm_lift_pos(const mm_lift_rec_t *rec, int64_t c)
{
	int z = mm_lift_zone(rec, c);
	if (!(z & 1)) return mm_lift_block_off(rec, z / 2) + c;               /* a block */
	return mm_lift_block_off(rec, z / 2) + mm_lift_gap_at(rec, z / 2).jct;   /* novel bases: the junction */
}

/* 1 iff contig offset c is at or inside the novel bases of a gap and the
 * alignment (ending at c_end) has no reference coverage beyond them */
static int at_tail_gap(const mm_lift_rec_t *rec, int64_t c, int64_t c_end)
{
	int32_t g, n = mm_lift_n_gaps(rec);
	for (g = 0; g < n; ++g) {
		mm_lift_gap_t x = mm_lift_gap_at(rec, g);
		int64_t e = (int64_t)x.jct + x.ins_len;
		if (x.jct <= c && c <= e && c_end <= e) return 1;
	}
	return 0;
}

/* CIGAR buffer; cb_push() merges equal adjacent ops and drops empty ones */
typedef struct {
	uint32_t *a;
	int32_t n, cap;
	int over;      /* capacity exceeded; n keeps counting the required size */
} cig_buf_t;

static void cb_push(cig_buf_t *b, int32_t len, int op)
{
	if (len <= 0) return;
	if (b->n > 0 && !b->over && CIG_OP(b->a[b->n - 1]) == (uint32_t)op) {
		b->a[b->n - 1] = CIG_MK(CIG_LEN(b->a[b->n - 1]) + len, op);
		return;
	}
	if (b->n > 0 && b->over) {
		/* cannot inspect the buffer any more; count conservatively */
		++b->n;
		return;
	}
	if (b->n >= b->cap) { b->over = 1; ++b->n; return; }
	b->a[b->n++] = CIG_MK(len, op);
}

int mm_lift_cigar(const mm_lift_rec_t *rec, int32_t ctg_len,
                  int64_t rs, const uint32_t *cigar, int32_t n_cigar,
                  int32_t qs, int32_t qe, int rev,
                  int64_t *rs_out, int32_t *qs_out, int32_t *qe_out,
                  uint32_t *cigar_out, int32_t cap, int32_t *n_out)
{
	cig_buf_t b;
	int64_t c = rs;                 /* current contig offset */
	int32_t i;
	int32_t lead_clip = 0, tail_clip = 0;   /* query bases turned into clipping */
	int started = 0;                /* have we emitted a ref-consuming op? */
	int64_t first_ref = 0;          /* base-ref coord of the first such op */
	int32_t cur_blk = -1;           /* block of the last anchored base */
	int64_t c_end;                  /* contig end of the alignment */
	int64_t ref_span = 0;           /* reference span emitted */
	(void)ctg_len;

	b.a = cigar_out; b.n = 0; b.cap = cap; b.over = 0;

	/* pass 1: reject an alignment with no reference span or lying wholly
	 * inside the novel bases of one gap */
	{
		int64_t cc = rs;
		for (i = 0; i < n_cigar; ++i)
			if (op_consumes_ref(CIG_OP(cigar[i]))) cc += CIG_LEN(cigar[i]);
		c_end = cc;
		if (c_end <= rs) return -2;
		{
			int z0 = mm_lift_zone(rec, rs);
			if ((z0 & 1) && c_end <= mm_lift_zone_end(rec, z0)) return -2;
		}
	}

	/* pass 2: rewrite */
	for (i = 0; i < n_cigar; ++i) {
		int op = CIG_OP(cigar[i]);
		int32_t l = CIG_LEN(cigar[i]);

		if (!op_consumes_ref(op)) {
			/* pass through; a leading I, or one at a gap with nothing
			 * anchored after it, becomes clipping */
			if (op == OP_I && !started) lead_clip += l;
			else if (op == OP_I && at_tail_gap(rec, c, c_end)) tail_clip += l;
			else cb_push(&b, l, op);
			continue;
		}

		/* split [c, c+l) across the zones */
		while (l > 0) {
			int z = mm_lift_zone(rec, c);
			int64_t zone_end = mm_lift_zone_end(rec, z);
			int32_t take;

			take = (int32_t)(zone_end - c < l ? zone_end - c : l);
			if (take <= 0) break;   /* cannot happen */

			if (!(z & 1)) {         /* block z/2 */
				int32_t blk = z / 2, g;
				if (started && blk > cur_blk) {
					/* deletions of the gaps crossed, covered on both sides */
					for (g = cur_blk; g < blk; ++g) {
						int32_t dl = mm_lift_gap_at(rec, g).del_len;
						if (dl > 0) { cb_push(&b, dl, OP_D); ref_span += dl; }
					}
				}
				if (!started) {
					started = 1;
					first_ref = mm_lift_block_off(rec, blk) + c;
				}
				cur_blk = blk;
				cb_push(&b, take, op);
				ref_span += take;
			} else {                /* novel bases of gap z/2 */
				if (op == OP_M || op == OP_EQ || op == OP_X) {
					/* aligned to novel sequence: insertion, or clipping at an end */
					if (!started)                 lead_clip += take;
					else if (c_end <= zone_end)   tail_clip += take;
					else                          cb_push(&b, take, OP_I);
				}
				/* D/N inside the gap: nothing to emit */
			}
			c += take;
			l -= take;
		}
	}

	/* trim leading/trailing D/N */
	if (!b.over) {
		int32_t lo = 0, hi = b.n;
		while (lo < hi && (CIG_OP(b.a[lo]) == OP_D || CIG_OP(b.a[lo]) == OP_N)) {
			first_ref += CIG_LEN(b.a[lo]);
			ref_span  -= CIG_LEN(b.a[lo]);
			++lo;
		}
		while (hi > lo && (CIG_OP(b.a[hi - 1]) == OP_D || CIG_OP(b.a[hi - 1]) == OP_N)) {
			ref_span -= CIG_LEN(b.a[hi - 1]);
			--hi;
		}
		/* a leading/trailing I that survived trimming becomes clipping */
		while (lo < hi && CIG_OP(b.a[lo]) == OP_I) { lead_clip += CIG_LEN(b.a[lo]); ++lo; }
		while (hi > lo && CIG_OP(b.a[hi - 1]) == OP_I) { tail_clip += CIG_LEN(b.a[hi - 1]); --hi; }

		if (lo > 0 && hi > lo) memmove(b.a, b.a + lo, (size_t)(hi - lo) * sizeof(uint32_t));
		b.n = hi - lo;
		if (b.n == 0) return -2;   /* nothing anchored survived */
	}

	*n_out = b.n;
	if (b.over) return -1;

	*rs_out = first_ref;
	/* clips are in reference order, qs/qe in forward-read order */
	if (!rev) {
		*qs_out = qs + lead_clip;
		*qe_out = qe - tail_clip;
	} else {
		*qs_out = qs + tail_clip;
		*qe_out = qe - lead_clip;
	}
	(void)ref_span;
	return 0;
}

/* same-locus test: shared diagonal on a common query interval */

typedef struct { int32_t qs, qe; int64_t diag; } diag_seg_t;

/* Match blocks of an alignment as (forward-read interval, diagonal); the
 * diagonal is t-q on the forward strand and t+q on the reverse. Returns the
 * number written, at most cap. */
static int32_t diag_segs(const mm_reg1_t *r, diag_seg_t *out, int32_t cap)
{
	int32_t n = 0;
	if (r->p == NULL || r->p->n_cigar == 0) {
		if (cap > 0) { out[0].qs = r->qs; out[0].qe = r->qe; out[0].diag = r->rev? (int64_t)r->rs + r->qe : (int64_t)r->rs - r->qs; n = 1; }
		return n;
	} else {
		int64_t t = r->rs; int32_t q = r->rev? r->qe : r->qs; uint32_t k;
		for (k = 0; k < r->p->n_cigar; ++k) {
			int op = CIG_OP(r->p->cigar[k]);
			int32_t l = CIG_LEN(r->p->cigar[k]);
			if (op == OP_M || op == OP_EQ || op == OP_X) {
				if (n < cap) {
					if (r->rev) { out[n].qs = q - l; out[n].qe = q; out[n].diag = t + q; }
					else        { out[n].qs = q; out[n].qe = q + l; out[n].diag = t - q; }
				}
				++n; t += l; q += r->rev? -l : l;
			} else if (op == OP_I || op == OP_S) {
				q += r->rev? -l : l;
			} else if (op == OP_D || op == OP_N) {
				t += l;
			}
		}
		return n < cap ? n : cap;
	}
}

int mm_lift_same_locus(const mm_reg1_t *a, const mm_reg1_t *b)
{
	diag_seg_t da[64], db[64];
	int32_t na, nb, i, j;
	if (a->rid != b->rid || a->rev != b->rev) return 0;
	na = diag_segs(a, da, 64);
	nb = diag_segs(b, db, 64);
	for (i = 0; i < na; ++i)
		for (j = 0; j < nb; ++j)
			if (da[i].diag == db[j].diag && da[i].qs < db[j].qe && db[j].qs < da[i].qe)
				return 1;
	return 0;
}
