/* --hap-dp: haplotype-consistent DP score (see mm_hapdp.h). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "kalloc.h"
#include "mmpriv.h"
#include "mm_hapdp.h"

static const int8_t hd_nt4[256] = {
	['A'] = 1, ['C'] = 2, ['G'] = 3, ['T'] = 4, ['a'] = 1, ['c'] = 2, ['g'] = 3, ['t'] = 4
}; /* code + 1; 0 = not ACGT */

/* run statistics, printed at exit (only when a table was loaded) */
static long hd_cand = 0, hd_cand_alt = 0, hd_walked = 0, hd_changed = 0, hd_strings = 0, hd_trunc = 0, hd_mid = 0, hd_mid_flip = 0, hd_base_mismatch = 0;
static int hd_loaded = 0;

__attribute__((destructor)) static void mm_hapdp_report(void)
{
	if (hd_loaded)
		fprintf(stderr, "[hap-dp] candidates with panel SNVs in the window %ld, with a window ALT in the read %ld, changed by a haplotype string %ld; strings %ld (scored %ld); "
				"windows truncated by --hap-dp-max %ld; middles below the ungapped bar with SNVs %ld, of them kept ungapped by a string %ld; "
				"base re-score != mm_update_extra %ld\n", hd_cand, hd_cand_alt, hd_changed, hd_strings, hd_walked, hd_trunc, hd_mid, hd_mid_flip, hd_base_mismatch);
}

void mm_hapdp_destroy(mm_hapdp_t *hd)
{
	if (hd == NULL) return;
	free(hd->s);
	free(hd->pool);
	free(hd);
}

static int hd_rec_cmp(const void *a_, const void *b_)
{
	const mm_hd_rec_t *a = (const mm_hd_rec_t*)a_, *b = (const mm_hd_rec_t*)b_;
	return a->pos != b->pos ? (a->pos < b->pos ? -1 : 1) : a->alt != b->alt ? (a->alt < b->alt ? -1 : 1) : 0;
}

typedef struct { int32_t rid; mm_hd_rec_t r; } hd_tmp_t;

static int hd_tmp_cmp(const void *a_, const void *b_)
{
	const hd_tmp_t *a = (const hd_tmp_t*)a_, *b = (const hd_tmp_t*)b_;
	if (a->rid != b->rid) return a->rid < b->rid ? -1 : 1;
	return hd_rec_cmp(&a->r, &b->r);
}

static int hd_parse_hex(const char *p, char **end, uint64_t *x)
{
	uint64_t v = 0;
	int n = 0;
	for (; *p; ++p, ++n) {
		int c = *p, d;
		if (c >= '0' && c <= '9') d = c - '0';
		else if (c >= 'a' && c <= 'f') d = c - 'a' + 10;
		else if (c >= 'A' && c <= 'F') d = c - 'A' + 10;
		else break;
		if (n >= 16) return -1;
		v = v << 4 | (uint64_t)d;
	}
	*end = (char*)p, *x = v;
	return n > 0 ? 0 : -1;
}

mm_hapdp_t *mm_hapdp_load(const mm_idx_t *mi, const char *fn)
{
	gzFile fp;
	mm_hapdp_t *hd;
	char *line, last_name[1024];
	int line_cap = 1 << 16, last_rid = -1, max_bit = -1;
	hd_tmp_t *t = 0;
	int64_t n_t = 0, m_t = 0, i, k;
	uint32_t *cnt;

	if (mi->h == 0) mm_idx_index_name((mm_idx_t*)mi);
	fp = strcmp(fn, "-") == 0 ? gzdopen(0, "r") : gzopen(fn, "r");
	if (fp == 0) return NULL;
	hd = (mm_hapdp_t*)calloc(1, sizeof(mm_hapdp_t));
	hd->n_seq = mi->n_seq;
	hd->max_str = 16;
	hd->n_hap = -1;
	line = (char*)malloc(line_cap);
	last_name[0] = 0;
	while (gzgets(fp, line, line_cap) != NULL) {
		char *p = line, *name, *end;
		long pos;
		int rid, l, c;
		uint64_t hi, lo;
		uint8_t b;
		l = strlen(line);
		if (l == line_cap - 1 && line[l-1] != '\n') {
			++hd->n_bad;
			while (gzgets(fp, line, line_cap) != NULL) { l = strlen(line); if (l < line_cap - 1 || line[l-1] == '\n') break; }
			continue;
		}
		while (l > 0 && (line[l-1] == '\n' || line[l-1] == '\r')) line[--l] = 0;
		if (l == 0) continue;
		if (line[0] == '#') {
			if (strncmp(line, "#n_hap\t", 7) == 0) hd->n_hap = atoi(line + 7);
			continue;
		}
		++hd->n_lines;
		name = p;
		while (*p && *p != '\t') ++p;
		if (*p != '\t') { ++hd->n_bad; continue; }
		*p++ = 0;
		pos = strtol(p, &end, 10);
		if (end == p || *end != '\t' || pos < 1) { ++hd->n_bad; continue; }
		p = end + 1;
		c = hd_nt4[(uint8_t)*p];
		if (c == 0 || p[1] != '\t') { ++hd->n_bad; continue; }
		p += 2;
		if (hd_parse_hex(p, &end, &hi) < 0 || *end != '\t') { ++hd->n_bad; continue; }
		p = end + 1;
		if (hd_parse_hex(p, &end, &lo) < 0 || (*end != 0 && *end != '\t')) { ++hd->n_bad; continue; }
		if (last_rid >= 0 && strcmp(name, last_name) == 0) rid = last_rid;
		else {
			rid = mm_idx_name2id(mi, name);
			if (strlen(name) < sizeof(last_name)) { strcpy(last_name, name); last_rid = rid; }
			else last_name[0] = 0, last_rid = -1;
		}
		if (rid < 0) { ++hd->n_unknown_seq; continue; }
		if ((uint64_t)pos > mi->seq[rid].len) { ++hd->n_out_of_range; continue; }
		mm_idx_getseq(mi, rid, pos - 1, pos, &b);
		if (b > 3 || b == c - 1) { ++hd->n_eq_ref; continue; }   /* ALT equals the base, or base is N */
		if (hi == 0 && lo == 0) { ++hd->n_nocarrier; continue; }
		if (n_t == m_t) {
			m_t = m_t ? m_t + (m_t >> 1) : 1 << 16;
			t = (hd_tmp_t*)realloc(t, m_t * sizeof(hd_tmp_t));
		}
		t[n_t].rid = rid;
		t[n_t].r.pos = (uint32_t)(pos - 1);
		t[n_t].r.alt = (uint32_t)(c - 1);
		t[n_t].r.m[0] = lo, t[n_t].r.m[1] = hi;
		++n_t;
		for (k = 0; k < 64; ++k) {
			if (lo >> k & 1) max_bit = max_bit > k ? max_bit : (int)k;
			if (hi >> k & 1) max_bit = max_bit > 64 + k ? max_bit : (int)(64 + k);
		}
	}
	gzclose(fp);
	free(line);
	if (n_t > 0) qsort(t, n_t, sizeof(hd_tmp_t), hd_tmp_cmp);
	/* merge equal (rid, pos, alt); a haplotype marked at two ALTs of one position keeps the first */
	for (i = 1, k = 0; i < n_t; ++i) {
		if (t[i].rid == t[k].rid && t[i].r.pos == t[k].r.pos && t[i].r.alt == t[k].r.alt) {
			t[k].r.m[0] |= t[i].r.m[0], t[k].r.m[1] |= t[i].r.m[1], ++hd->n_dup;
		} else t[++k] = t[i];
	}
	n_t = n_t ? k + 1 : 0;
	for (i = 0; i < n_t; ++i) {
		int64_t j;
		for (j = i - 1; j >= 0 && t[j].rid == t[i].rid && t[j].r.pos == t[i].r.pos; --j) {
			uint64_t o0 = t[i].r.m[0] & t[j].r.m[0], o1 = t[i].r.m[1] & t[j].r.m[1];
			if (o0 || o1) { ++hd->n_conflict; t[i].r.m[0] &= ~o0, t[i].r.m[1] &= ~o1; }
		}
	}
	for (i = 0, k = 0; i < n_t; ++i)   /* a record emptied by the conflict rule goes */
		if (t[i].r.m[0] || t[i].r.m[1]) t[k++] = t[i];
	n_t = k;
	if (hd->n_hap < 0) hd->n_hap = max_bit + 1;
	else if (max_bit >= hd->n_hap) {
		fprintf(stderr, "[ERROR] --hap-dp: a mask sets haplotype %d but #n_hap is %d\n", max_bit, hd->n_hap);
		free(t); mm_hapdp_destroy(hd);
		return NULL;
	}
	hd->n_rec = n_t;
	hd->pool = (mm_hd_rec_t*)malloc((n_t ? n_t : 1) * sizeof(mm_hd_rec_t));
	hd->s = (mm_hd_seq_t*)calloc(mi->n_seq, sizeof(mm_hd_seq_t));
	cnt = (uint32_t*)calloc(mi->n_seq, sizeof(uint32_t));
	for (i = 0; i < n_t; ++i) hd->pool[i] = t[i].r, ++cnt[t[i].rid];
	for (i = 0, k = 0; i < (int64_t)mi->n_seq; ++i) {
		hd->s[i].n = cnt[i];
		hd->s[i].a = hd->pool + k;
		k += cnt[i];
	}
	free(cnt);
	free(t);
	hd_loaded = 1;
	return hd;
}

int mm_hapdp_load_mask(mm_hapdp_t *hd, const char *fn)
{
	FILE *fp;
	char buf[4096];
	int n = 0;
	if ((fp = fopen(fn, "r")) == NULL) return -1;
	hd->allow[0] = hd->allow[1] = 0;
	while (fgets(buf, sizeof(buf), fp)) {
		char *p = buf, *end;
		if (buf[0] == '#') continue;
		for (;;) {
			long h;
			while (*p == ' ' || *p == '\t' || *p == '\n' || *p == '\r' || *p == ',') ++p;
			if (*p == 0) break;
			h = strtol(p, &end, 10);
			if (end == p || h < 0 || h >= 128 || (hd->n_hap > 0 && h >= hd->n_hap)) { fclose(fp); return -2; }
			hd->allow[h >> 6] |= 1ULL << (h & 63);
			++n, p = end;
		}
	}
	fclose(fp);
	hd->has_allow = 1;
	return n > 0 ? 0 : -3;
}

static inline uint32_t hd_lower_bound(const mm_hd_seq_t *s, uint32_t x)
{
	uint32_t lo = 0, hi = s->n;
	while (lo < hi) {
		uint32_t mid = lo + ((hi - lo) >> 1);
		if (s->a[mid].pos < x) lo = mid + 1;
		else hi = mid;
	}
	return lo;
}

int32_t mm_hapdp_range(const mm_hapdp_t *hd, int32_t rid, int32_t st, int32_t en, uint32_t *first)
{
	const mm_hd_seq_t *s;
	uint32_t j, e;
	if (hd == NULL || rid < 0 || (uint32_t)rid >= hd->n_seq || st >= en || st < 0) return 0;
	s = &hd->s[rid];
	if (s->n == 0) return 0;
	j = hd_lower_bound(s, (uint32_t)st);
	if (j == s->n || s->a[j].pos >= (uint32_t)en) return 0;
	for (e = j; e < s->n && s->a[e].pos < (uint32_t)en; ++e) {}
	*first = j;
	return (int32_t)(e - j);
}

void mm_hapdp_win_free(void *km, mm_hd_win_t *w)
{
	if (w->on) kfree(km, w->on);
	if (w->cnt) kfree(km, w->cnt);
	w->on = 0, w->cnt = 0;
}

int32_t mm_hapdp_window(void *km, const mm_hapdp_t *hd, int32_t rid, int32_t st, int32_t en, mm_hd_win_t *w)
{
	uint32_t first;
	int32_t n_rec, j, s, n_str, cap, a, b;
	uint64_t u[2], (*cl)[2];
	memset(w, 0, sizeof(*w));
	n_rec = mm_hapdp_range(hd, rid, st, en, &first);
	if (n_rec == 0) return 0;
	w->n_rec = n_rec;
	w->rec = hd->s[rid].a + first;
	u[0] = u[1] = 0;
	for (j = 0; j < n_rec; ++j) u[0] |= w->rec[j].m[0], u[1] |= w->rec[j].m[1];
	if (hd->has_allow) u[0] &= hd->allow[0], u[1] &= hd->allow[1];
	cap = __builtin_popcountll(u[0]) + __builtin_popcountll(u[1]);
	w->n_hap_win = cap;
	if (cap == 0) return 0;
	/* distinct strings: partition the carriers u by every record's mask */
	cl = (uint64_t(*)[2])kmalloc(km, (size_t)cap * 16);
	cl[0][0] = u[0], cl[0][1] = u[1], n_str = 1;
	for (j = 0; j < n_rec && n_str < cap; ++j) {
		const uint64_t *m = w->rec[j].m;
		int32_t n0 = n_str;
		for (s = 0; s < n0; ++s) {
			uint64_t x0 = cl[s][0] & m[0], x1 = cl[s][1] & m[1], y0 = cl[s][0] & ~m[0], y1 = cl[s][1] & ~m[1];
			if ((x0 | x1) && (y0 | y1)) {
				cl[s][0] = x0, cl[s][1] = x1;
				cl[n_str][0] = y0, cl[n_str][1] = y1, ++n_str;
			}
		}
	}
	w->on = (uint8_t*)kmalloc(km, (size_t)n_str * n_rec);
	w->cnt = (int32_t*)kmalloc(km, (size_t)n_str * sizeof(int32_t));
	{	/* order by count, descending; ties by the lowest haplotype */
		int32_t *ord = (int32_t*)kmalloc(km, (size_t)n_str * 2 * sizeof(int32_t)), *lo = ord + n_str;
		for (s = 0; s < n_str; ++s) {
			ord[s] = s;
			lo[s] = cl[s][0] ? __builtin_ctzll(cl[s][0]) : 64 + __builtin_ctzll(cl[s][1]);
		}
		for (a = 1; a < n_str; ++a) {
			int32_t x = ord[a], cx = __builtin_popcountll(cl[x][0]) + __builtin_popcountll(cl[x][1]);
			for (b = a - 1; b >= 0; --b) {
				int32_t y = ord[b], cy = __builtin_popcountll(cl[y][0]) + __builtin_popcountll(cl[y][1]);
				if (cy > cx || (cy == cx && lo[y] < lo[x])) break;
				ord[b + 1] = y;
			}
			ord[b + 1] = x;
		}
		for (s = 0; s < n_str; ++s) {
			int32_t x = ord[s], h = lo[x];
			w->cnt[s] = __builtin_popcountll(cl[x][0]) + __builtin_popcountll(cl[x][1]);
			for (j = 0; j < n_rec; ++j) w->on[(size_t)s * n_rec + j] = w->rec[j].m[h >> 6] >> (h & 63) & 1;
		}
		kfree(km, ord);
	}
	kfree(km, cl);
	if (n_str > hd->max_str - 1) {   /* keep the most frequent strings; the base takes one slot */
		int32_t keep = hd->max_str - 1 > 0 ? hd->max_str - 1 : 0;
		w->truncated = n_str - keep;
		n_str = keep;
	}
	w->n_str = n_str;
	return n_str;
}

static inline int32_t hd_sc(const mm_mapopt_t *opt, int cq, int ct)
{
	if (cq > 3 || ct > 3) return opt->sc_ambi > 0 ? -opt->sc_ambi : opt->sc_ambi;
	return cq == ct ? opt->a : -opt->b;
}

int32_t mm_hapdp_ungapped_best(void *km, const mm_hapdp_t *hd, const mm_mapopt_t *opt, int32_t rid, int32_t st, int32_t len,
							   const uint8_t *qseq, const uint8_t *tseq, int32_t base_score)
{
	mm_hd_win_t w;
	uint32_t first;
	int32_t s, j, best = base_score, n_rec, any = 0;
	const mm_hd_rec_t *rec;
	n_rec = mm_hapdp_range(hd, rid, st, st + len, &first);
	if (n_rec == 0) return base_score;
	__atomic_fetch_add(&hd_mid, 1, __ATOMIC_RELAXED);
	rec = hd->s[rid].a + first;
	for (j = 0; j < n_rec && !any; ++j) any = (qseq[rec[j].pos - st] == rec[j].alt);
	if (!any) return base_score;   /* a string can only gain where the read carries its ALT */
	if (mm_hapdp_window(km, hd, rid, st, st + len, &w) == 0) { mm_hapdp_win_free(km, &w); return base_score; }
	for (s = 0; s < w.n_str; ++s) {   /* ungapped: every position contributes on its own */
		const uint8_t *on = w.on + (size_t)s * w.n_rec;
		int32_t sc = base_score;
		for (j = 0; j < w.n_rec; ++j) {
			int32_t o;
			if (!on[j]) continue;
			o = (int32_t)w.rec[j].pos - st;
			sc += hd_sc(opt, qseq[o], w.rec[j].alt) - hd_sc(opt, qseq[o], tseq[o]);
		}
		if (sc > best) best = sc;
	}
	mm_hapdp_win_free(km, &w);
	return best;
}

/* Score the target with substitutions sub[] (offset<<2|base, ascending) along the CIGAR as mm_update_extra()
 * does. Returns dp_max; *sum = sum of substitution scores over M/=/X, *mlen = matches, *n_alt = substituted
 * positions matched by the read. */
static int32_t hd_walk(const mm_reg1_t *r, const int8_t *mat, int8_t q, int8_t e, int log_gap, const uint8_t *qseq, const uint8_t *tseq,
					   int32_t n_sub, const uint32_t *sub, int64_t *sum, int32_t *mlen, int32_t *n_alt)
{
	const mm_extra_t *p = r->p;
	uint32_t k, l;
	int32_t qoff = 0, toff = 0, si = 0, ml = 0, na = 0;
	int64_t sm = 0;
	double s = 0.0, max = 0.0;
	for (k = 0; k < p->n_cigar; ++k) {
		uint32_t op = p->cigar[k]&0xf, len = p->cigar[k]>>4;
		if (op == MM_CIGAR_MATCH || op == MM_CIGAR_EQ_MATCH || op == MM_CIGAR_X_MISMATCH) {
			for (l = 0; l < len; ++l) {
				int cq = qseq[qoff + l], ct = tseq[toff + l], t = toff + l;
				while (si < n_sub && (int32_t)(sub[si] >> 2) < t) ++si;
				if (si < n_sub && (int32_t)(sub[si] >> 2) == t) {
					ct = sub[si] & 3;
					if (cq == ct) ++na;
				}
				if (ct <= 3 && cq <= 3 && ct == cq) ++ml;
				sm += mat[ct * 5 + cq];
				s += mat[ct * 5 + cq];
				if (s < 0) s = 0;
				else max = max > s ? max : s;
			}
			qoff += len, toff += len;
		} else if (op == MM_CIGAR_INS || op == MM_CIGAR_DEL) {
			if (log_gap) s -= q + (double)e * mg_log2(1.0 + len);
			else s -= q + e;
			if (s < 0) s = 0;
			if (op == MM_CIGAR_INS) qoff += len; else toff += len;
		} else if (op == MM_CIGAR_N_SKIP) toff += len;
	}
	*sum = sm, *mlen = ml, *n_alt = na;
	return (int32_t)(max + .499);
}

int mm_hapdp_rescore(void *km, const mm_idx_t *mi, const mm_mapopt_t *opt, const int8_t *mat, int log_gap,
					 const uint8_t *qseq, mm_reg1_t *r)
{
	const mm_hapdp_t *hd = mi->hapdp;
	mm_hd_win_t w;
	uint8_t *tseq, *qb;
	uint32_t *sub, first, k;
	int32_t s, j, n_sub, best_s = -1, best_max, best_ml, best_na = 0, ml0, na0, mx, ml, na, n_rec, any = 0, qoff = 0, toff = 0;
	int64_t sum0, sum, best_sum;
	const mm_hd_rec_t *rec;
	r->hd = 0, r->hd_dmlen = 0;
	if (hd == NULL || r->p == NULL || r->p->n_cigar == 0 || r->inv || r->re <= r->rs) return 0;
	n_rec = mm_hapdp_range(hd, r->rid, r->rs, r->re, &first);
	if (n_rec == 0) return 0;
	__atomic_fetch_add(&hd_cand, 1, __ATOMIC_RELAXED);
	/* read base aligned to each record (4: deleted/skipped); a string with none of its ALTs in the read cannot
	 * beat the base and is skipped */
	rec = hd->s[r->rid].a + first;
	qb = (uint8_t*)kmalloc(km, n_rec);
	memset(qb, 4, n_rec);
	for (k = 0, j = 0; k < r->p->n_cigar && j < n_rec; ++k) {
		uint32_t op = r->p->cigar[k]&0xf, len = r->p->cigar[k]>>4;
		if (op == MM_CIGAR_MATCH || op == MM_CIGAR_EQ_MATCH || op == MM_CIGAR_X_MISMATCH) {
			for (; j < n_rec && rec[j].pos < (uint32_t)(r->rs + toff + len); ++j)
				if (rec[j].pos >= (uint32_t)(r->rs + toff)) qb[j] = qseq[qoff + (rec[j].pos - r->rs - toff)];
			qoff += len, toff += len;
		} else if (op == MM_CIGAR_INS) qoff += len;
		else if (op == MM_CIGAR_DEL || op == MM_CIGAR_N_SKIP) {
			for (; j < n_rec && rec[j].pos < (uint32_t)(r->rs + toff + len); ++j) {}
			toff += len;
		}
	}
	for (j = 0; j < n_rec && !any; ++j) any = (qb[j] == rec[j].alt);
	if (!any) { kfree(km, qb); return 0; }
	__atomic_fetch_add(&hd_cand_alt, 1, __ATOMIC_RELAXED);
	if (mm_hapdp_window(km, hd, r->rid, r->rs, r->re, &w) == 0) { mm_hapdp_win_free(km, &w); kfree(km, qb); return 0; }
	__atomic_fetch_add(&hd_strings, w.n_str, __ATOMIC_RELAXED);
	if (w.truncated) __atomic_fetch_add(&hd_trunc, 1, __ATOMIC_RELAXED);
	tseq = (uint8_t*)kmalloc(km, r->re - r->rs);
	sub = (uint32_t*)kmalloc(km, (w.n_rec + 1) * sizeof(uint32_t));
	mm_idx_getseq(mi, r->rid, r->rs, r->re, tseq);
	best_max = hd_walk(r, mat, opt->q, opt->e, log_gap, qseq, tseq, 0, sub, &sum0, &ml0, &na0);
	if (best_max != r->p->dp_max) __atomic_fetch_add(&hd_base_mismatch, 1, __ATOMIC_RELAXED);
	best_sum = sum0, best_ml = ml0;
	for (s = 0; s < w.n_str; ++s) {
		const uint8_t *on = w.on + (size_t)s * w.n_rec;
		int32_t gain = 0;
		for (j = 0; j < w.n_rec; ++j) gain += on[j] && qb[j] == w.rec[j].alt;
		if (gain == 0) continue;
		__atomic_fetch_add(&hd_walked, 1, __ATOMIC_RELAXED);
		for (j = 0, n_sub = 0; j < w.n_rec; ++j)
			if (on[j]) sub[n_sub++] = (w.rec[j].pos - (uint32_t)r->rs) << 2 | w.rec[j].alt;
		mx = hd_walk(r, mat, opt->q, opt->e, log_gap, qseq, tseq, n_sub, sub, &sum, &ml, &na);
		if (mx > best_max || (mx == best_max && sum > best_sum))
			best_s = s, best_max = mx, best_sum = sum, best_ml = ml, best_na = na;
	}
	kfree(km, sub);
	kfree(km, tseq);
	kfree(km, qb);
	mm_hapdp_win_free(km, &w);
	if (best_s < 0) return 0;   /* the base is best or tied */
	r->p->dp_max = best_max;
	r->p->dp_score += (int32_t)(best_sum - sum0);
	r->hd_dmlen = best_ml - ml0;
	r->mlen += r->hd_dmlen;
	r->hd = best_na;
	__atomic_fetch_add(&hd_changed, 1, __ATOMIC_RELAXED);
	return 1;
}

void mm_hapdp_count_mid_flip(void) { __atomic_fetch_add(&hd_mid_flip, 1, __ATOMIC_RELAXED); }
