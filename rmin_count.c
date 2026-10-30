/* rmin_count: count read minimizers (mm_sketch, no HPC) and canonical k-mers for rmin.
 *   rmin_count [-k 27] [-w 14] [-n MAXREADS_PER_FILE] [-M] -o PREFIX reads.fq.gz ...
 * Output PREFIX.mm.rmin / PREFIX.km.rmin (RMIN1, little-endian): "RMIN1\0\0\0", uint32 k, uint32 w
 * (1 for k-mers), uint64 n_distinct, n_occurrences, n_reads; then n x uint64 key<<8|min(count,255),
 * ascending, key = hash64(canonical k-mer). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/resource.h>
#include "bseq.h"
#include "minimap.h"
#include "mmpriv.h"
#include "kvec.h"
#include "ketopt.h"

static inline uint64_t hash64(uint64_t key, uint64_t mask) /* as in sketch.c */
{
	key = (~key + (key << 21)) & mask;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)) & mask;
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)) & mask;
	key = key ^ key >> 28;
	key = (key + (key << 31)) & mask;
	return key;
}

typedef struct { size_t n, m; uint64_t *a; } u64v_t;

static void push(u64v_t *v, uint64_t x)
{
	if (v->n == v->m) {
		v->m = v->m? v->m + (v->m >> 1) : 1<<20;
		v->a = (uint64_t*)realloc(v->a, v->m * 8);
		if (v->a == 0) { fprintf(stderr, "[E::rmin_count] out of memory at %zu keys\n", v->n); exit(1); }
	}
	v->a[v->n++] = x;
}

static void all_kmers(const char *s, int len, int k, u64v_t *v)
{
	uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, f = 0, r = 0;
	int i, l = 0;
	for (i = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)s[i]];
		if (c < 4) {
			f = (f << 2 | c) & mask;
			r = (r >> 2) | (3ULL ^ c) << shift1;
			if (++l >= k && f != r) push(v, hash64(f < r? f : r, mask));
		} else l = 0;
	}
}

static int write_set(const char *fn, u64v_t *v, int k, int w, uint64_t n_reads)
{
	size_t i, j, n_out = 0;
	uint64_t n_occ = v->n;
	FILE *fp;
	radix_sort_64(v->a, v->a + v->n);
	for (i = 0, j = 0; i < v->n; i = j) { /* compact in place: key<<8 | min(cnt,255) */
		uint64_t cnt;
		for (j = i + 1; j < v->n && v->a[j] == v->a[i]; ++j);
		cnt = j - i;
		v->a[n_out++] = v->a[i] << 8 | (cnt > 255? 255 : cnt);
	}
	if ((fp = fopen(fn, "wb")) == 0) { fprintf(stderr, "[E::rmin_count] cannot write %s\n", fn); return -1; }
	{
		uint32_t kw[2] = { (uint32_t)k, (uint32_t)w };
		uint64_t h[3] = { n_out, n_occ, n_reads };
		fwrite("RMIN1\0\0\0", 1, 8, fp);
		fwrite(kw, 4, 2, fp);
		fwrite(h, 8, 3, fp);
		if (fwrite(v->a, 8, n_out, fp) != n_out) { fprintf(stderr, "[E::rmin_count] write error %s\n", fn); fclose(fp); return -1; }
	}
	fclose(fp);
	fprintf(stderr, "[M::rmin_count] %s: k=%d w=%d occurrences=%lu distinct=%zu\n", fn, k, w, (unsigned long)n_occ, n_out);
	return 0;
}

int main(int argc, char *argv[])
{
	int c, k = 27, w = 14, no_kmers = 0, f;
	int64_t max_reads = -1;
	const char *prefix = 0;
	uint64_t n_reads = 0, n_bases = 0;
	u64v_t vm = {0,0,0}, vk = {0,0,0};
	mm128_v mv = {0,0,0};
	ketopt_t o = KETOPT_INIT;
	double t0 = (double)time(0);
	char fn[4096];
	struct rusage ru;

	while ((c = ketopt(&o, argc, argv, 1, "k:w:n:Mo:", 0)) >= 0) {
		if (c == 'k') k = atoi(o.arg);
		else if (c == 'w') w = atoi(o.arg);
		else if (c == 'n') max_reads = atoll(o.arg);
		else if (c == 'M') no_kmers = 1;
		else if (c == 'o') prefix = o.arg;
	}
	if (prefix == 0 || o.ind >= argc || k < 2 || k > 28 || w < 1 || w > 255) {
		fprintf(stderr, "Usage: rmin_count [-k 27] [-w 14] [-n maxReadsPerFile] [-M] -o PREFIX in1.fq.gz [in2.fq.gz ...]\n");
		return 1;
	}
	for (f = o.ind; f < argc; ++f) {
		mm_bseq_file_t *fp = mm_bseq_open(argv[f]);
		int64_t n_this = 0;
		if (fp == 0) { fprintf(stderr, "[E::rmin_count] cannot open %s\n", argv[f]); return 1; }
		for (;;) {
			int n = 0, i;
			mm_bseq1_t *seq = mm_bseq_read(fp, 200000000, 0, &n);
			if (seq == 0 || n == 0) { free(seq); break; }
			for (i = 0; i < n; ++i) {
				if (max_reads < 0 || n_this < max_reads) {
					mv.n = 0;
					if (seq[i].l_seq > 0) mm_sketch(0, seq[i].seq, seq[i].l_seq, w, k, 0, 0, &mv);
					for (size_t j = 0; j < mv.n; ++j) push(&vm, mv.a[j].x >> 8);
					if (!no_kmers) all_kmers(seq[i].seq, seq[i].l_seq, k, &vk);
					++n_this, n_bases += seq[i].l_seq;
				}
				free(seq[i].seq); free(seq[i].name);
			}
			free(seq);
			if (max_reads >= 0 && n_this >= max_reads) break;
		}
		mm_bseq_close(fp);
		n_reads += n_this;
		fprintf(stderr, "[M::rmin_count] %s: reads %ld (total %lu, bases %lu), minimizer occ %zu, kmer occ %zu, %.0f s\n",
				argv[f], (long)n_this, (unsigned long)n_reads, (unsigned long)n_bases, vm.n, vk.n, (double)time(0) - t0);
	}
	free(mv.a);
	snprintf(fn, sizeof(fn), "%s.mm.rmin", prefix);
	if (write_set(fn, &vm, k, w, n_reads) != 0) return 1;
	free(vm.a);
	if (!no_kmers) {
		snprintf(fn, sizeof(fn), "%s.km.rmin", prefix);
		if (write_set(fn, &vk, k, 1, n_reads) != 0) return 1;
		free(vk.a);
	}
	getrusage(RUSAGE_SELF, &ru);
	fprintf(stderr, "[M::rmin_count] done: %.0f s wall, peak RSS %.2f GB\n", (double)time(0) - t0, ru.ru_maxrss / 1048576.0);
	return 0;
}
