/* Index-side half of the projection: building the per-sequence table and
 * applying it to a candidate set.
 *
 * Kept separate from mm_lift.c so that the projection arithmetic stays free of
 * any dependency on the rest of minimap2 and test/test_mm_lift.c can link
 * against mm_lift.c alone.  Everything here needs mm_idx_t / mm_reg1_t.
 */

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
	free(lt);
}

int mm_lift_build(mm_idx_t *mi, const char *fn_map, int32_t jct, int64_t *n_unresolved)
{
	mm_lift_t *lt;
	size_t n;
	int32_t i;
	int64_t n_ok = 0, n_bad = 0;
	int32_t *delta = NULL, *ins3 = NULL, *del3 = NULL;

	if (mi == NULL || mi->n_seq == 0 || jct <= 0) return -1;
	n = (size_t)mi->n_seq;

	/* reuse minimap2's own name->rid table (index.c:143,162) rather than
	 * building a second one */
	if (mm_idx_index_name(mi) != 0) return -1;

	lt = (mm_lift_t*)calloc(1, sizeof(mm_lift_t));
	if (lt == NULL) return -1;
	lt->n = (int32_t)n;
	lt->a = (mm_lift_rec_t*)calloc(n, sizeof(mm_lift_rec_t));
	delta = (int32_t*)calloc(n, sizeof(int32_t));
	ins3  = (int32_t*)calloc(n, sizeof(int32_t));
	del3  = (int32_t*)calloc(n, sizeof(int32_t));
	if (lt->a == NULL || delta == NULL || ins3 == NULL || del3 == NULL) { free(delta); free(ins3); free(del3); free(lt->a); free(lt); return -1; }
	for (i = 0; i < lt->n; ++i) delta[i] = INT32_MIN;   /* "not supplied" */

	if (fn_map != NULL) {
		FILE *fp = fopen(fn_map, "r");
		char buf[1024];
		if (fp == NULL) {
			fprintf(stderr, "[ERROR] --lift-map: cannot open %s\n", fn_map);
			free(delta); free(ins3); free(del3); free(lt->a); free(lt);
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
				delta[id] = (int32_t)v1;                 /* two columns: d = |ALT|-|REF| */
			} else {                                     /* three columns: ins_len, del_len */
				v2 = strtol(tab2 + 1, NULL, 10);
				ins3[id] = (int32_t)v1;
				del3[id] = (int32_t)v2;
				delta[id] = INT32_MIN + 1;               /* marker: explicit pair supplied */
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

		/* a contig closer than jct-1 to the chromosome start has a shortened
		 * left flank: the generator takes lf = min(flank, pos-1), so the
		 * junction is min(jct, pos) and off_a is then 0 */
		jct_i = pos < jct ? (int32_t)pos : jct;
		d = delta[i];
		if (d == INT32_MIN + 1) {                        /* explicit (ins_len, del_len) */
			if (ins3[i] < 0 || del3[i] < 0 || (ins3[i] == 0 && del3[i] == 0)) { ++n_bad; continue; }
			lt->a[i].base_rid = base_id;
			lt->a[i].jct      = jct_i;
			lt->a[i].off_a    = pos - jct_i;
			lt->a[i].ins_len  = ins3[i];
			lt->a[i].del_len  = del3[i];
			++n_ok;
			continue;
		}
		if (d == INT32_MIN) {
			/* Not supplied.  len = 200 + |ALT|, and an insertion has |REF| = 1,
			 * so d = len - 201.  The generator requires |ALT| >= 5, hence a
			 * genuine insertion has d >= 4.  Anything shorter is a deletion
			 * whose |REF| the length does not carry: leave it unresolved
			 * rather than guess.  See the metadata audit in mm_lift.h. */
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

	free(delta); free(ins3); free(del3);
	mm_lift_destroy((mm_lift_t*)mi->lift);
	mi->lift = (struct mm_lift_s*)lt;
	if (n_unresolved) *n_unresolved = n_bad;
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
		if ((int32_t)r->p->n_cigar + 4 > cap) {
			cap = (int32_t)r->p->n_cigar + 4;
			buf = (uint32_t*)malloc((size_t)cap * 4);
			if (buf == NULL) continue;
		}

		rc = mm_lift_cigar(rec, (int32_t)mi->seq[r->rid].len,
		                   r->rs, r->p->cigar, (int32_t)r->p->n_cigar,
		                   r->qs, r->qe, r->rev,
		                   &rs_out, &qs_out, &qe_out, buf, cap, &n_out);
		if (rc == 0) {
			/* grow r->p if needed; same accounting as align.c:313-316 --
			 * capacity counts uint32 words and includes the header */
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
		/* rc == -1 cannot happen: cap >= n_cigar+4 and the splice adds at
		 * most 2 ops.  rc == -2 means no base-reference representation
		 * exists; leave the reg exactly as it was. */
		if (buf != stack_buf) free(buf);
	}
	return n_done;
}
