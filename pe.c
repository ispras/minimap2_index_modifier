#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mmpriv.h"
#include "kvec.h"

void mm_select_sub_multi(void *km, float pri_ratio, float pri1, float pri2, int max_gap_ref, int min_diff, int best_n, int n_segs, const int *qlens, int *n_, mm_reg1_t *r)
{
	if (pri_ratio > 0.0f && *n_ > 0) {
		int i, k, n = *n_, n_2nd = 0;
		int max_dist = n_segs == 2? qlens[0] + qlens[1] + max_gap_ref : 0;
		for (i = k = 0; i < n; ++i) {
			int to_keep = 0;
			if (r[i].parent == i) { // primary
				to_keep = 1;
			} else if (r[i].on_ctg) {
				/* --lift: keep contig chains; a fragment chain on the base would
				 * outscore them before DP compares the alignments */
				to_keep = 2;
			} else if (r[i].score + min_diff >= r[r[i].parent].score) {
				to_keep = 1;
			} else {
				mm_reg1_t *p = &r[r[i].parent], *q = &r[i];
				if (p->rev == q->rev && p->rid == q->rid && q->re - p->rs < max_dist && p->re - q->rs < max_dist) { // child and parent are close on the ref
					if (q->score >= p->score * pri1)
						to_keep = 1;
				} else {
					int is_par_both = (n_segs == 2 && p->qs < qlens[0] && p->qe > qlens[0]);
					int is_chi_both = (n_segs == 2 && q->qs < qlens[0] && q->qe > qlens[0]);
					if (is_chi_both || is_chi_both == is_par_both) {
						if (q->score >= p->score * pri_ratio)
							to_keep = 1;
					} else { // the remaining case: is_chi_both == 0 && is_par_both == 1
						if (q->score >= p->score * pri2)
							to_keep = 1;
					}
				}
			}
			if (to_keep == 1 && r[i].parent != i) {
				if (n_2nd++ >= best_n) to_keep = 0; // don't keep if there are too many secondary hits
			}
			if (to_keep) r[k++] = r[i];
			else if (r[i].p) free(r[i].p);
		}
		if (k != n) mm_sync_regs(km, k, r); // removing hits requires sync()
		*n_ = k;
	}
}

void mm_set_pe_thru(const int *qlens, int *n_regs, mm_reg1_t **regs)
{
	int s, i, n_pri[2], pri[2];
	n_pri[0] = n_pri[1] = 0;
	pri[0] = pri[1] = -1;
	for (s = 0; s < 2; ++s)
		for (i = 0; i < n_regs[s]; ++i)
			if (regs[s][i].id == regs[s][i].parent)
				++n_pri[s], pri[s] = i;
	if (n_pri[0] == 1 && n_pri[1] == 1) {
		mm_reg1_t *p = &regs[0][pri[0]];
		mm_reg1_t *q = &regs[1][pri[1]];
		if (p->rid == q->rid && p->rev == q->rev && abs(p->rs - q->rs) < 3 && abs(p->re - q->re) < 3
			&& ((p->qs == 0 && qlens[1] - q->qe == 0) || (q->qs == 0 && qlens[0] - p->qe == 0)))
		{
			p->pe_thru = q->pe_thru = 1;
		}
	}
}

#include "ksort.h"

typedef struct {
	int s, rev;
	uint64_t key;
	mm_reg1_t *r;
} pair_arr_t;

#define sort_key_pair(a) ((a).key)
KRADIX_SORT_INIT(pair, pair_arr_t, sort_key_pair, 8)

/* Log base lambda of the scoring matrix (as vg's recover_log_base()): the root of
 * 1/4 e^{lambda a} + 3/4 e^{-lambda b} = 1, assuming GC 0.5. */
double mm_pair_log_base(int a, int b)
{
	static int ca = -1, cb = -1;
	static double cl = 0.0;
	double lo = 1e-6, hi = 10.0, lam, f;
	int it;
	if (a == ca && b == cb) return cl;
	for (it = 0; it < 200; ++it) {
		lam = 0.5 * (lo + hi);
		f = 0.25 * exp(lam * a) + 0.75 * exp(-lam * b) - 1.0;
		if (f < 0.0) lo = lam; else hi = lam;
	}
	cl = 0.5 * (lo + hi), cb = b, ca = a;   /* benign race: every thread computes the same value */
	return cl;
}

/* fragment length of a candidate pair, as in mm_mate_resolve(): reference span of the two ends
 * + clipped bases at the two outer ends + (insertions - deletions) of both */
static int pe_frag_on_read(const mm_reg1_t *x, int sx, const mm_reg1_t *y, int sy, const int *qlens)
{
	const mm_reg1_t *left = x->rs <= y->rs ? x : y, *right = x->rs <= y->rs ? y : x;
	int ql_left = qlens[left == x ? sx : sy], ql_right = qlens[right == x ? sx : sy];
	int lo = x->rs < y->rs ? x->rs : y->rs, hi = x->re > y->re ? x->re : y->re, sp = hi - lo, t;
	sp += left->rev ? ql_left - left->qe : left->qs;            /* clip at the reference-left end */
	sp += right->rev ? right->qs : ql_right - right->qe;        /* clip at the reference-right end */
	for (t = 0; t < 2; ++t) {
		const mm_reg1_t *c = t == 0 ? x : y;
		if (c->p) { uint32_t k; for (k = 0; k < c->p->n_cigar; ++k) { int op = c->p->cigar[k]&0xf, l = c->p->cigar[k]>>4; if (op == MM_CIGAR_INS) sp += l; else if (op == MM_CIGAR_DEL) sp -= l; } }
	}
	return sp;
}

static long pmq_pairs = 0, pmq_multi = 0, pmq_prior_moved = 0, pmq_lowered = 0, pmq_raised = 0;
__attribute__((destructor)) static void pmq_report(void)
{
	if (pmq_pairs > 0)
		fprintf(stderr, "[pair-mapq] fragments with a pair %ld (more than one pair %ld); best pair changed by the fragment-length prior %ld; end mapq lowered %ld, raised %ld\n",
				pmq_pairs, pmq_multi, pmq_prior_moved, pmq_lowered, pmq_raised);
}

void mm_pair(void *km, int max_gap_ref, int pe_bonus, int sub_diff, int match_sc, const int *qlens, int *n_regs, mm_reg1_t **regs, int no_blend)
{
	mm_pair_core(km, max_gap_ref, pe_bonus, sub_diff, match_sc, qlens, n_regs, regs, no_blend, 0, 0.0, 0.0, 0.0);
}

/* pair_mode 0: default pairing.  pair_mode 1 (vg) / 2 (vgup), as vg giraffe:
 *  - pair score ps = dp_max0 + dp_max1 - (d - mu)^2 / (2 sd^2 lambda), not below the worse end's dp_max;
 *    no prior when mu or sd is 0; the pair with the largest ps is chosen;
 *  - raw = -10 log10(1 - e^{lambda ps_best} / sum_i e^{lambda ps_i}) over pairs within pe_bonus of the
 *    best sum; a single pair gets a null alignment of score 0;
 *  - both ends get mapq (int)(min(raw, 120) / 2); vgup: an end takes it only if higher than its own. */
void mm_pair_core(void *km, int max_gap_ref, int pe_bonus, int sub_diff, int match_sc, const int *qlens, int *n_regs, mm_reg1_t **regs, int no_blend,
				  int pair_mode, double frag_mu, double frag_sd, double lambda)
{
	int i, j, s, n, last[2], dp_thres, segs = 0, max_idx[2], raw_idx[2] = {-1, -1};
	int64_t max, raw_max = -1;
	double ps_max = 0.0;
	int best_k = -1;
	pair_arr_t *a;
	kvec_t(uint64_t) sc = {0,0,0};
	kvec_t(double) psv = {0,0,0};

	a = (pair_arr_t*)kmalloc(km, (n_regs[0] + n_regs[1]) * sizeof(pair_arr_t));
	for (s = n = 0, dp_thres = 0; s < 2; ++s) {
		int max = 0;
		for (i = 0; i < n_regs[s]; ++i) {
			a[n].s = s;
			a[n].r = &regs[s][i];
			a[n].rev = a[n].r->rev;
			a[n].key = (uint64_t)a[n].r->rid << 32 | a[n].r->rs<<1 | (s^a[n].rev);
			max = max > a[n].r->p->dp_max? max : a[n].r->p->dp_max;
			++n;
			segs |= 1<<s;
		}
		dp_thres += max;
	}
	if (segs != 3) {
		kfree(km, a); // only one end is mapped
		return;
	}
	dp_thres -= pe_bonus;
	if (dp_thres < 0) dp_thres = 0;
	radix_sort_pair(a, a + n);

	max = -1;
	max_idx[0] = max_idx[1] = -1;
	last[0] = last[1] = -1;
	kv_resize(uint64_t, km, sc, (size_t)n);
	for (i = 0; i < n; ++i) {
		if (a[i].key & 1) { // reverse first read or forward second read
			mm_reg1_t *q, *r;
			if (last[a[i].rev] < 0) continue;
			r = a[i].r;
			q = a[last[a[i].rev]].r;
			if (r->rid != q->rid || r->rs - q->re > max_gap_ref) continue;
			for (j = last[a[i].rev]; j >= 0; --j) {
				int64_t score;
				if (a[j].rev != a[i].rev || a[j].s == a[i].s) continue;
				q = a[j].r;
				if (r->rid != q->rid || r->rs - q->re > max_gap_ref) break;
				if (r->p->dp_max + q->p->dp_max < dp_thres) continue;
				score = (int64_t)(r->p->dp_max + q->p->dp_max) << 32 | (r->hash + q->hash);
				if (pair_mode == 0) {
					if (score > max)
						max = score, max_idx[a[j].s] = j, max_idx[a[i].s] = i;
				} else {
					double ps = (double)(r->p->dp_max + q->p->dp_max), worse = (double)(r->p->dp_max < q->p->dp_max ? r->p->dp_max : q->p->dp_max);
					if (frag_mu > 0.0 && frag_sd > 0.0) {
						double dev = (double)pe_frag_on_read(q, a[j].s, r, a[i].s, qlens) - frag_mu;
						ps -= dev * dev / (2.0 * frag_sd * frag_sd) / lambda;
						if (ps < worse) ps = worse;
					}
					if (max < 0 || ps > ps_max || (ps == ps_max && score > max))
						max = score, ps_max = ps, best_k = (int)psv.n, max_idx[a[j].s] = j, max_idx[a[i].s] = i;
					if (score > raw_max) raw_max = score, raw_idx[a[j].s] = j, raw_idx[a[i].s] = i;
					kv_push(double, km, psv, ps);
				}
				kv_push(uint64_t, km, sc, score);
			}
		} else { // forward first read or reverse second read
			last[a[i].rev] = i;
		}
	}
	if (sc.n > 1)
		radix_sort_64(sc.a, sc.a + sc.n);

	if (sc.n > 0 && max > 0) { // found at least one pair
		int n_sub = 0, mapq_pe;
		mm_reg1_t *r[2];
		r[0] = a[max_idx[0]].r, r[1] = a[max_idx[1]].r;
		r[0]->proper_frag = r[1]->proper_frag = 1;
		for (s = 0; s < 2; ++s) {
			if (r[s]->id != r[s]->parent) { // then lift to primary and update parent
				mm_reg1_t *p = &regs[s][r[s]->parent];
				for (i = 0; i < n_regs[s]; ++i)
					if (regs[s][i].parent == p->id)
						regs[s][i].parent = r[s]->id;
				p->mapq = 0;
			}
			if (!r[s]->sam_pri) { // then sync sam_pri
				for (i = 0; i < n_regs[s]; ++i)
					regs[s][i].sam_pri = 0;
				r[s]->sam_pri = 1;
			}
		}
		if (pair_mode != 0) {
			double T = 0.0, raw;
			int q, rawi;
			for (i = 0; i < (int)psv.n; ++i)
				if (i != best_k) T += exp(lambda * (psv.a[i] - ps_max));
			if (psv.n == 1) T = exp(-lambda * ps_max);   /* one pair: vg adds a null alignment of score 0 */
			raw = T > 0.0 ? 10.0 * log10(1.0 + 1.0 / T) : 2147483647.0;
			rawi = raw >= 2147483647.0 ? 2147483647 : (int)raw;
			q = (int)((rawi < 120 ? rawi : 120) / 2.0);
			__atomic_fetch_add(&pmq_pairs, 1, __ATOMIC_RELAXED);
			if (psv.n > 1) __atomic_fetch_add(&pmq_multi, 1, __ATOMIC_RELAXED);
			if (raw_idx[0] != max_idx[0] || raw_idx[1] != max_idx[1]) __atomic_fetch_add(&pmq_prior_moved, 1, __ATOMIC_RELAXED);
			for (s = 0; s < 2; ++s) {
				if (pair_mode == 2 && r[s]->mapq >= q) continue;
				if (r[s]->mapq > q) __atomic_fetch_add(&pmq_lowered, 1, __ATOMIC_RELAXED);
				else if (r[s]->mapq < q) __atomic_fetch_add(&pmq_raised, 1, __ATOMIC_RELAXED);
				r[s]->mapq = q;
			}
			goto pair_done;
		}
		mapq_pe = r[0]->mapq > r[1]->mapq? r[0]->mapq : r[1]->mapq;
		for (i = 0; i < (int)sc.n; ++i)
			if ((sc.a[i]>>32) + sub_diff >= (uint64_t)max>>32)
				++n_sub;
		if (sc.n > 1) {
			int mapq_pe_alt;
			mapq_pe_alt = (int)(6.02f * ((max>>32) - (sc.a[sc.n - 2]>>32)) / match_sc - 4.343f * logf(n_sub)); // n_sub > 0 because it counts the optimal, too
			mapq_pe = mapq_pe < mapq_pe_alt? mapq_pe : mapq_pe_alt;
		}
		if (no_blend) { /* --mate-resolve: take the pair-level mapq, no 0.2/0.8 blend */
			if (r[0]->mapq < mapq_pe) r[0]->mapq = mapq_pe;
			if (r[1]->mapq < mapq_pe) r[1]->mapq = mapq_pe;
		} else {
			if (r[0]->mapq < mapq_pe) r[0]->mapq = (int)(.2f * r[0]->mapq + .8f * mapq_pe + .499f);
			if (r[1]->mapq < mapq_pe) r[1]->mapq = (int)(.2f * r[1]->mapq + .8f * mapq_pe + .499f);
		}
		if (sc.n == 1) {
			if (r[0]->mapq < 2) r[0]->mapq = 2;
			if (r[1]->mapq < 2) r[1]->mapq = 2;
		} else if ((uint64_t)max>>32 > sc.a[sc.n - 2]>>32) {
			if (r[0]->mapq < 1) r[0]->mapq = 1;
			if (r[1]->mapq < 1) r[1]->mapq = 1;
		}
	}
pair_done:

	kfree(km, a);
	kfree(km, sc.a);
	if (psv.a) kfree(km, psv.a);

	mm_set_pe_thru(qlens, n_regs, regs);
}
