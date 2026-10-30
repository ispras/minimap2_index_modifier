/* Building the projection table and applying it to mapped regs. The pure
 * arithmetic lives in mm_lift.c, which has no other minimap2 dependency. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "mmpriv.h"
#include "mm_lift.h"

void mm_lift_destroy(mm_lift_t *lt)
{
	if (lt == NULL) return;
	free(lt->a);
	free(lt->xg);
	free(lt);
}

/* Parse "jct:ins:del;jct:ins:del..." into pool[*n_pool..], growing it.
 * Returns the number of gaps parsed, or -1 if malformed or too many. */
static int parse_xgaps(const char *s, mm_lift_gap_t **pool, int64_t *n_pool, int64_t *m_pool)
{
	int k = 0;
	for (;;) {
		long v[3]; int f; char *end;
		for (f = 0; f < 3; ++f) {
			if (*s < '0' || *s > '9') return -1;
			v[f] = strtol(s, &end, 10);
			s = end;
			if (f < 2) { if (*s != ':') return -1; ++s; }
		}
		if (k + 1 >= MM_LIFT_MAX_GAPS) return -1;   /* plus gap 0 */
		if (*n_pool == *m_pool) {
			int64_t m = *m_pool ? *m_pool * 2 : 1024;
			mm_lift_gap_t *np = (mm_lift_gap_t*)realloc(*pool, (size_t)m * sizeof(mm_lift_gap_t));
			if (np == NULL) return -1;
			*pool = np; *m_pool = m;
		}
		(*pool)[*n_pool].jct = (int32_t)v[0];
		(*pool)[*n_pool].ins_len = (int32_t)v[1];
		(*pool)[*n_pool].del_len = (int32_t)v[2];
		++*n_pool; ++k;
		if (*s == ';') { ++s; continue; }
		if (*s == '\n' || *s == '\r' || *s == '\0') return k;
		return -1;
	}
}

int mm_lift_build(mm_idx_t *mi, const char *fn_map, int32_t jct, int64_t *n_unresolved)
{
	mm_lift_t *lt;
	size_t n;
	int32_t i;
	int64_t n_ok = 0, n_bad = 0;
	int64_t n_jct4 = 0;   /* contigs whose jct came from the map file */
	int32_t *delta = NULL, *ins3 = NULL, *del3 = NULL, *jct4 = NULL;
	/* extra gaps per rid: pool[xoff, xoff+xn); xn = -1: malformed */
	int64_t *xoff = NULL; int32_t *xn = NULL;
	mm_lift_gap_t *pool = NULL; int64_t n_pool = 0, m_pool = 0, n_multi = 0;

	if (mi == NULL || mi->n_seq == 0 || jct <= 0) return -1;
	n = (size_t)mi->n_seq;

	if (mm_idx_index_name(mi) != 0) return -1;

	lt = (mm_lift_t*)calloc(1, sizeof(mm_lift_t));
	if (lt == NULL) return -1;
	lt->n = (int32_t)n;
	lt->a = (mm_lift_rec_t*)calloc(n, sizeof(mm_lift_rec_t));
	delta = (int32_t*)calloc(n, sizeof(int32_t));
	ins3  = (int32_t*)calloc(n, sizeof(int32_t));
	del3  = (int32_t*)calloc(n, sizeof(int32_t));
	jct4  = (int32_t*)calloc(n, sizeof(int32_t));   /* 0 = not supplied */
	xoff  = (int64_t*)calloc(n, sizeof(int64_t));
	xn    = (int32_t*)calloc(n, sizeof(int32_t));   /* 0 = no extra gaps */
	if (lt->a == NULL || delta == NULL || ins3 == NULL || del3 == NULL || jct4 == NULL || xoff == NULL || xn == NULL) { free(delta); free(ins3); free(del3); free(jct4); free(xoff); free(xn); free(lt->a); free(lt); return -1; }
	for (i = 0; i < lt->n; ++i) delta[i] = INT32_MIN;   /* "not supplied" */

	if (fn_map != NULL) {
		FILE *fp = fopen(fn_map, "r");
		char buf[1024];
		if (fp == NULL) {
			fprintf(stderr, "[ERROR] --lift-map: cannot open %s\n", fn_map);
			free(delta); free(ins3); free(del3); free(jct4); free(xoff); free(xn); free(lt->a); free(lt);
			return -1;
		}
		while (fgets(buf, sizeof buf, fp)) {
			char *tab = strchr(buf, '\t'), *tab2, *end;
			int id;
			long v1, v2;
			if (tab == NULL) continue;
			*tab = '\0';
			id = mm_idx_name2id(mi, buf);
			if (id < 0 || id >= lt->n) continue;
			v1 = strtol(tab + 1, &end, 10);
			tab2 = (*end == '\t') ? end : NULL;
			if (tab2 == NULL) {
				delta[id] = (int32_t)v1;                 /* d = |ALT|-|REF| */
			} else {                                     /* ins_len, del_len [, jct [, extra gaps]] */
				char *end2;
				v2 = strtol(tab2 + 1, &end2, 10);
				ins3[id] = (int32_t)v1;
				del3[id] = (int32_t)v2;
				xn[id] = 0;                              /* a repeated line replaces the earlier one */
				delta[id] = INT32_MIN + 1;               /* marker: explicit pair supplied */
				if (*end2 == '\t') {                     /* jct */
					char *end3;
					long v3 = strtol(end2 + 1, &end3, 10);
					if (v3 > 0) { jct4[id] = (int32_t)v3; ++n_jct4; }
					if (*end3 == '\t') {                 /* extra gaps */
						int64_t o = n_pool;
						int k = parse_xgaps(end3 + 1, &pool, &n_pool, &m_pool);
						if (k < 0) { n_pool = o; xn[id] = -1; }
						else { xoff[id] = o; xn[id] = k; }
					}
				}
			}
		}
		fclose(fp);
	}

	for (i = 0; i < lt->n; ++i) {
		int32_t chrom_len, d, jct_i;
		int64_t pos;
		int base_id;
		char save, *nm = mi->seq[i].name;

		lt->a[i].base_rid = -1;
		if (nm == NULL) continue;
		if (!mm_lift_parse_name(nm, &chrom_len, &pos)) continue;   /* a base sequence */

		/* resolve the chromosome prefix against the index's own name table */
		save = nm[chrom_len];
		nm[chrom_len] = '\0';
		base_id = mm_idx_name2id(mi, nm);
		nm[chrom_len] = save;
		if (base_id < 0) { ++n_bad; continue; }

		/* junction from the map file, else min(jct, pos) (short left flank at a
		 * chromosome start) */
		if (jct4[i] > 0) {
			jct_i = jct4[i];
			if ((int64_t)jct_i > pos || jct_i >= (int32_t)mi->seq[i].len) { ++n_bad; continue; }   /* impossible geometry */
		} else jct_i = pos < jct ? (int32_t)pos : jct;
		d = delta[i];
		if (d == INT32_MIN + 1) {                        /* explicit (ins_len, del_len) */
			if (ins3[i] < 0 || del3[i] < 0 || (ins3[i] == 0 && del3[i] == 0)) { ++n_bad; continue; }
			if (jct4[i] > 0 && jct_i + ins3[i] >= (int32_t)mi->seq[i].len) { ++n_bad; continue; }   /* no block B left */
			if (xn[i] < 0) { ++n_bad; continue; }   /* malformed extra gaps */
			if (xn[i] > 0) {
				/* requires an explicit jct; gaps ascend with non-empty blocks
				 * between and after them */
				int32_t k, prev_end = jct_i + ins3[i], ok = (jct4[i] > 0);
				for (k = 0; ok && k < xn[i]; ++k) {
					const mm_lift_gap_t *x = &pool[xoff[i] + k];
					if (x->ins_len < 0 || x->del_len < 0 || (x->ins_len == 0 && x->del_len == 0)) ok = 0;
					else if (x->jct <= prev_end) ok = 0;
					else prev_end = x->jct + x->ins_len;
				}
				if (!ok || prev_end >= (int32_t)mi->seq[i].len) { ++n_bad; continue; }
				lt->a[i].n_xgap = xn[i];   /* xgap pointer set once the pool is final */
				++n_multi;
			}
			lt->a[i].base_rid = base_id;
			lt->a[i].jct      = jct_i;
			lt->a[i].off_a    = pos - jct_i;
			lt->a[i].ins_len  = ins3[i];
			lt->a[i].del_len  = del3[i];
			++n_ok;
			continue;
		}
		if (xn[i] != 0) { ++n_bad; continue; }   /* extra gaps need the explicit form */
		if (d == INT32_MIN) {
			/* not supplied: infer an insertion d = len - 201 (100 bp flanks,
			 * |REF| = 1, |ALT| >= 5); a deletion cannot be inferred */
			int32_t inferred = (int32_t)mi->seq[i].len - 201;
			if (inferred >= 4) d = inferred;
			else { ++n_bad; continue; }
		}
		if (d == 0) { ++n_bad; continue; }   /* neither insertion nor deletion */

		lt->a[i].base_rid = base_id;
		lt->a[i].jct      = jct_i;
		lt->a[i].off_a    = pos - jct_i;
		lt->a[i].ins_len  = d > 0 ?  d : 0;
		lt->a[i].del_len  = d < 0 ? -d : 0;
		++n_ok;
	}

	/* the pool is final: point the haplotype contigs into it */
	lt->xg = pool; lt->n_xg = n_pool; lt->n_multi = n_multi;
	for (i = 0; i < lt->n; ++i)
		if (lt->a[i].base_rid >= 0 && lt->a[i].n_xgap > 0) lt->a[i].xgap = pool + xoff[i];
		else lt->a[i].n_xgap = 0, lt->a[i].xgap = NULL;

	free(delta); free(ins3); free(del3); free(jct4); free(xoff); free(xn);
	mm_lift_destroy((mm_lift_t*)mi->lift);
	mi->lift = (struct mm_lift_s*)lt;
	if (n_unresolved) *n_unresolved = n_bad;
	lt->n_jct_from_map = n_jct4;
	return (int)n_ok;
}

int mm_lift_regs(const mm_idx_t *mi, int n_regs, mm_reg1_t *regs)
{
	const mm_lift_t *lt = (const mm_lift_t*)mi->lift;
	uint32_t stack_buf[128];
	int i, n_done = 0;

	if (lt == NULL || regs == NULL) return 0;

	for (i = 0; i < n_regs; ++i) {
		mm_reg1_t *r = &regs[i];
		const mm_lift_rec_t *rec;
		int64_t rs_out;
		int32_t qs_out, qe_out, n_out, cap;
		uint32_t *buf;
		int rc;

		if (r->rid < 0 || r->rid >= lt->n) continue;
		rec = &lt->a[r->rid];
		if (rec->base_rid < 0) continue;

		if (r->p == NULL || r->p->n_cigar == 0) {
			/* chaining only, no CIGAR: project the endpoints */
			int64_t a = mm_lift_pos(rec, r->rs);
			int64_t b = mm_lift_pos(rec, r->re);
			if (b < a) b = a;
			r->lift_rid = r->rid; r->lift_rs = r->rs;
			r->rid = rec->base_rid;
			r->rs = (int32_t)a;
			r->re = (int32_t)b;
			r->lifted = 1;
			++n_done;
			continue;
		}

		cap = (int32_t)(sizeof(stack_buf) / sizeof(stack_buf[0]));
		buf = stack_buf;
		/* the splice adds at most three ops per gap */
		if ((int32_t)r->p->n_cigar + 3 * mm_lift_n_gaps(rec) + 1 > cap) {
			cap = (int32_t)r->p->n_cigar + 3 * mm_lift_n_gaps(rec) + 1;
			buf = (uint32_t*)malloc((size_t)cap * 4);
			if (buf == NULL) continue;
		}

		rc = mm_lift_cigar(rec, (int32_t)mi->seq[r->rid].len,
		                   r->rs, r->p->cigar, (int32_t)r->p->n_cigar,
		                   r->qs, r->qe, r->rev,
		                   &rs_out, &qs_out, &qe_out, buf, cap, &n_out);
		if (rc == 0) {
			/* grow r->p; capacity is in uint32 words including the header */
			if ((uint32_t)n_out + sizeof(mm_extra_t)/4 > r->p->capacity) {
				uint32_t ncap = (uint32_t)n_out + sizeof(mm_extra_t)/4;
				mm_extra_t *np;
				kroundup32(ncap);
				np = (mm_extra_t*)realloc(r->p, (size_t)ncap * 4);
				if (np == NULL) { if (buf != stack_buf) free(buf); continue; }
				r->p = np;
				r->p->capacity = ncap;
			}
			memcpy(r->p->cigar, buf, (size_t)n_out * 4);
			r->p->n_cigar = (uint32_t)n_out;
			r->lift_rid = r->rid; r->lift_rs = r->rs;
			r->rid = rec->base_rid;
			r->rs  = (int32_t)rs_out;
			r->re  = (int32_t)(rs_out + mm_lift_ref_span(buf, n_out));
			r->qs  = qs_out;
			r->qe  = qe_out;
			r->lifted = 1;
			++n_done;
		}
		/* on rc == -2 (no base representation) the reg is left unchanged */
		if (buf != stack_buf) free(buf);
	}
	return n_done;
}
