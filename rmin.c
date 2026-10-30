/* rmin: read-based filter of modifier-added minimizers (see rmin.h). Called
 * from add_variants() in sketch.c before each candidate is pushed. */
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "rmin.h"
#include "khash.h"

KHASH_SET_INIT_STR(rstr)

rmin_opt_t rmin_opt = { RMIN_NONE, 2, 0.5, 0.0, 11, 0, 0, 0, 0, 0 };

typedef struct { uint64_t *a; size_t n; uint32_t k, w; uint64_t n_occ, n_reads; } rmin_set_t;

static rmin_set_t S_mm, S_km, S_ref;
static khash_t(rstr) *S_list;
static FILE *fp_dump;
static uint64_t combo_id;

/* statistics; sequence class: 0 base, 1 contig, 2 haplotype contig (_h);
 * minimizer class: 0 overlaps a substituted ALT, 1 reference-only k-mer */
static uint64_t st_refcoll[3], st_cand[3][2], st_kept[3][2], st_combo[3], st_combo_c[3], st_combo_lost[3], st_combo_nocur[3],
	st_combo_partial[3], st_rec[3], st_rec_c[3], st_rec_lost[3], st_cnt_missing[3][2];
static int cur_cls = -1; static uint64_t rec_cand, rec_kept;

int rmin_parse_mode(const char *s)
{
	if (strcmp(s, "strict") == 0) return RMIN_STRICT;
	if (strcmp(s, "combo") == 0) return RMIN_COMBO;
	if (strcmp(s, "kmer") == 0) return RMIN_KMER;
	if (strcmp(s, "rand") == 0) return RMIN_RAND;
	if (strcmp(s, "list") == 0) return RMIN_LIST;
	return -1;
}

int rmin_active(void) { return rmin_opt.mode != RMIN_NONE; }

static int load_set(const char *fn, rmin_set_t *s)
{
	FILE *fp = fopen(fn, "rb");
	char magic[8];
	uint32_t kw[2];
	uint64_t h[3];
	if (fp == 0) { fprintf(stderr, "[E::rmin] cannot open %s: %s\n", fn, strerror(errno)); return -1; }
	if (fread(magic, 1, 8, fp) != 8 || memcmp(magic, "RMIN1\0\0\0", 8) != 0) { fprintf(stderr, "[E::rmin] %s: not an RMIN1 file\n", fn); fclose(fp); return -1; }
	if (fread(kw, 4, 2, fp) != 2 || fread(h, 8, 3, fp) != 3) { fprintf(stderr, "[E::rmin] %s: truncated header\n", fn); fclose(fp); return -1; }
	s->k = kw[0], s->w = kw[1], s->n = h[0], s->n_occ = h[1], s->n_reads = h[2];
	s->a = (uint64_t*)malloc(s->n * 8 + 8);
	if (fread(s->a, 8, s->n, fp) != s->n) { fprintf(stderr, "[E::rmin] %s: truncated body\n", fn); fclose(fp); return -1; }
	fclose(fp);
	for (size_t i = 1; i < s->n; ++i)
		if ((s->a[i] >> 8) <= (s->a[i-1] >> 8)) { fprintf(stderr, "[E::rmin] %s: keys not strictly sorted at %zu\n", fn, i); return -1; }
	fprintf(stderr, "[M::rmin] loaded %s: k=%u w=%u distinct=%lu occurrences=%lu reads=%lu\n", fn, s->k, s->w,
			(unsigned long)s->n, (unsigned long)s->n_occ, (unsigned long)s->n_reads);
	return 0;
}

static int set_count(const rmin_set_t *s, uint64_t key) /* -1: set not loaded; else count 0..255 */
{
	size_t lo = 0, hi;
	if (s->a == 0) return -1;
	hi = s->n;
	while (lo < hi) {
		size_t mid = lo + (hi - lo) / 2;
		if ((s->a[mid] >> 8) < key) lo = mid + 1; else hi = mid;
	}
	return (lo < s->n && (s->a[lo] >> 8) == key)? (int)(s->a[lo] & 0xff) : 0;
}

int rmin_init(void)
{
	rmin_opt_t *o = &rmin_opt;
	if (o->mode == RMIN_NONE) {
		if (o->fn_mm || o->fn_km || o->fn_list || o->fn_dump || o->fn_ref) { fprintf(stderr, "[E::rmin] --read-min-* files given without --read-min-mode\n"); return -1; }
		return 0;
	}
	if (o->fn_mm && load_set(o->fn_mm, &S_mm) != 0) return -1;
	if (o->fn_km && load_set(o->fn_km, &S_km) != 0) return -1;
	if (o->fn_ref && load_set(o->fn_ref, &S_ref) != 0) return -1;
	if (S_ref.a && S_ref.w != 1) { fprintf(stderr, "[E::rmin] --read-min-ref-kmers got a minimizer file (w=%u)\n", S_ref.w); return -1; }
	if ((o->mode == RMIN_STRICT || o->mode == RMIN_COMBO) && S_mm.a == 0) { fprintf(stderr, "[E::rmin] mode strict/combo needs --read-minimizers\n"); return -1; }
	if (o->mode == RMIN_KMER && S_km.a == 0) { fprintf(stderr, "[E::rmin] mode kmer needs --read-kmers\n"); return -1; }
	if (o->mode == RMIN_RAND && !(o->rate >= 0.0 && o->rate <= 1.0)) { fprintf(stderr, "[E::rmin] --read-min-rate must be in [0,1]\n"); return -1; }
	if (S_mm.a && S_mm.w == 1) { fprintf(stderr, "[E::rmin] --read-minimizers got a k-mer file (w=1)\n"); return -1; }
	if (S_km.a && S_km.w != 1) { fprintf(stderr, "[E::rmin] --read-kmers got a minimizer file (w=%u)\n", S_km.w); return -1; }
	if (o->mode == RMIN_LIST) {
		FILE *fp; char line[1024]; int absent; uint64_t n = 0;
		if (o->fn_list == 0 || (fp = fopen(o->fn_list, "r")) == 0) { fprintf(stderr, "[E::rmin] mode list needs a readable --read-min-list\n"); return -1; }
		S_list = kh_init(rstr);
		while (fgets(line, sizeof(line), fp)) {
			size_t l = strlen(line);
			while (l > 0 && (line[l-1] == '\n' || line[l-1] == '\r')) line[--l] = 0;
			if (l == 0 || line[0] == '#') continue;
			kh_put(rstr, S_list, strdup(line), &absent);
			++n;
		}
		fclose(fp);
		fprintf(stderr, "[M::rmin] loaded list %s: %lu lines, %u distinct\n", o->fn_list, (unsigned long)n, kh_size(S_list));
	}
	if (o->fn_dump) {
		if ((fp_dump = fopen(o->fn_dump, "w")) == 0) { fprintf(stderr, "[E::rmin] cannot write %s\n", o->fn_dump); return -1; }
		fprintf(fp_dump, "#combo\tchr\tcur\tsnvs\tend0\tstrand\thash\taltk\tmmcnt\tkmcnt\tkeep\n");
	}
	fprintf(stderr, "[M::rmin] mode=%d C=%d F=%g R=%g seed=%lu ref-collision-filter=%s\n", o->mode, o->min_cnt, o->frac, o->rate, (unsigned long)o->seed, o->fn_ref? o->fn_ref : "no");
	return 0;
}

static int chr_class(const char *chr)
{
	const char *p = strstr(chr, "_contig_");
	if (p == 0) return 0;
	return strstr(p, "_h")? 2 : 1;
}

void rmin_rec_begin(const char *chr, unsigned long curr_pos)
{
	if (!rmin_active()) return;
	cur_cls = chr_class(chr);
	rec_cand = rec_kept = 0;
	(void)curr_pos;
}

void rmin_rec_end(void)
{
	if (!rmin_active() || cur_cls < 0) return;
	++st_rec[cur_cls];
	if (rec_cand > 0) {
		++st_rec_c[cur_cls];
		if (rec_kept == 0) ++st_rec_lost[cur_cls];
	}
	cur_cls = -1;
}

static inline uint64_t splitmix64(uint64_t z)
{
	z += 0x9E3779B97F4A7C15ULL;
	z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
	z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
	return z ^ (z >> 31);
}

int rmin_filter(const mm128_t *a, int n, uint8_t *keep, const char *chr, unsigned long curr_pos,
                char **REF_arr, char **ALT_arr, unsigned long *POS_all, int N_SNP, int k)
{
	int i, j, n_kept = 0, has_cur = 0, cls = chr_class(chr);
	int mmc[64], kmc[64], altk[64];
	rmin_opt_t *o = &rmin_opt;
	(void)REF_arr;
	if (n > 64) { fprintf(stderr, "[E::rmin] %d candidates in one combination (>64)\n", n); exit(1); }
	++combo_id;
	++st_combo[cls];
	for (j = 0; j < N_SNP; ++j) if (POS_all[j] == curr_pos) has_cur = 1;
	if (!has_cur) ++st_combo_nocur[cls];
	for (i = 0; i < n; ++i) {
		uint64_t key = a[i].x >> 8;
		int64_t end0 = (int64_t)((uint32_t)a[i].y >> 1);
		mmc[i] = set_count(&S_mm, key);
		kmc[i] = set_count(&S_km, key);
		altk[i] = 0;
		for (j = 0; j < N_SNP; ++j) {
			int64_t p0 = (int64_t)POS_all[j] - 1;
			if (p0 <= end0 && p0 > end0 - k) { altk[i] = 1; break; }
		}
	}
	if (o->mode == RMIN_STRICT) {
		for (i = 0; i < n; ++i) keep[i] = (mmc[i] >= o->min_cnt);
	} else if (o->mode == RMIN_KMER) {
		for (i = 0; i < n; ++i) keep[i] = (kmc[i] >= o->min_cnt);
	} else if (o->mode == RMIN_COMBO) {
		int found = 0, all;
		for (i = 0; i < n; ++i) found += (mmc[i] >= o->min_cnt);
		all = (n == 0 || (double)found >= o->frac * n - 1e-9);
		for (i = 0; i < n; ++i) keep[i] = all;
	} else if (o->mode == RMIN_RAND) {
		for (i = 0; i < n; ++i) {
			uint64_t u = splitmix64(a[i].x ^ splitmix64(a[i].y ^ o->seed));
			keep[i] = ((double)(u >> 11) * (1.0 / 9007199254740992.0) >= o->rate);
		}
	} else if (o->mode == RMIN_LIST) {
		char buf[512]; int ok = 1;
		for (j = 0; j < N_SNP && ok; ++j) {
			snprintf(buf, sizeof(buf), "%s\t%lu\t%s", chr, POS_all[j], ALT_arr[j]);
			if (kh_get(rstr, S_list, buf) == kh_end(S_list)) ok = 0;
		}
		for (i = 0; i < n; ++i) keep[i] = ok;
	} else {
		for (i = 0; i < n; ++i) keep[i] = 1;
	}
	if (S_ref.a) /* drop ALT k-mers present in the reference */
		for (i = 0; i < n; ++i)
			if (altk[i] && set_count(&S_ref, a[i].x >> 8) > 0) { if (keep[i]) ++st_refcoll[cls]; keep[i] = 0; }
	for (i = 0; i < n; ++i) {
		int kc = altk[i]? 0 : 1;
		++st_cand[cls][kc];
		if (keep[i]) ++st_kept[cls][kc], ++n_kept;
		if (mmc[i] == 0) ++st_cnt_missing[cls][kc];
	}
	if (n > 0) {
		++st_combo_c[cls];
		if (n_kept == 0) ++st_combo_lost[cls];
		else if (n_kept < n) ++st_combo_partial[cls];
	}
	rec_cand += n, rec_kept += n_kept;
	if (fp_dump) {
		char snvs[2048]; int l = 0;
		snvs[0] = 0;
		for (j = 0; j < N_SNP && l < (int)sizeof(snvs) - 64; ++j)
			l += snprintf(snvs + l, sizeof(snvs) - l, "%s%lu:%s", j? "," : "", POS_all[j], ALT_arr[j]);
		for (i = 0; i < n; ++i)
			fprintf(fp_dump, "%lu\t%s\t%lu\t%s\t%u\t%u\t%lu\t%d\t%d\t%d\t%d\n", (unsigned long)combo_id, chr, curr_pos, snvs,
					(uint32_t)a[i].y >> 1, (unsigned)(a[i].y & 1), (unsigned long)(a[i].x >> 8), altk[i], mmc[i], kmc[i], keep[i]);
	}
	return n_kept;
}

void rmin_report(FILE *fp)
{
	static const char *cn[3] = { "base", "contig", "hapcontig" };
	int c;
	uint64_t tc = 0, tk = 0;
	if (!rmin_active()) return;
	for (c = 0; c < 3; ++c) {
		fprintf(fp, "[M::rmin_stat] %s: records %lu (with candidates %lu, lost all %lu); combos %lu (with candidates %lu, lost all %lu, partially kept %lu, without current SNV %lu); "
				"minimizers alt-kmer kept %lu of %lu (not in read minimizers %lu), ref-kmer kept %lu of %lu (not in read minimizers %lu)\n",
				cn[c], (unsigned long)st_rec[c], (unsigned long)st_rec_c[c], (unsigned long)st_rec_lost[c],
				(unsigned long)st_combo[c], (unsigned long)st_combo_c[c], (unsigned long)st_combo_lost[c], (unsigned long)st_combo_partial[c], (unsigned long)st_combo_nocur[c],
				(unsigned long)st_kept[c][0], (unsigned long)st_cand[c][0], (unsigned long)st_cnt_missing[c][0],
				(unsigned long)st_kept[c][1], (unsigned long)st_cand[c][1], (unsigned long)st_cnt_missing[c][1]);
		if (S_ref.a) fprintf(fp, "[M::rmin_stat] %s: alt-kmer minimizers dropped by reference collision (after mode) %lu\n", cn[c], (unsigned long)st_refcoll[c]);
		tc += st_cand[c][0] + st_cand[c][1], tk += st_kept[c][0] + st_kept[c][1];
	}
	fprintf(fp, "[M::rmin_stat] total: added minimizers kept %lu of %lu (dropped %lu, fraction dropped %.6f)\n",
			(unsigned long)tk, (unsigned long)tc, (unsigned long)(tc - tk), tc? (double)(tc - tk) / tc : 0.0);
	if (fp_dump) { fclose(fp_dump); fp_dump = 0; }
}

void rmin_destroy(void)
{
	free(S_mm.a); free(S_km.a); free(S_ref.a);
	S_mm.a = S_km.a = S_ref.a = 0;
	if (S_list) {
		khint_t it;
		for (it = kh_begin(S_list); it != kh_end(S_list); ++it)
			if (kh_exist(S_list, it)) free((char*)kh_key(S_list, it));
		kh_destroy(rstr, S_list);
		S_list = 0;
	}
}
