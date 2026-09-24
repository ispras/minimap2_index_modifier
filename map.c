#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <errno.h>
#include "kthread.h"
#include "kvec.h"
#include "kalloc.h"
#include "sdust.h"
#include "mmpriv.h"
#include "mm_lift.h"
#include <limits.h>
#include "bseq.h"
#include "khash.h"

mm_tbuf_t *mm_tbuf_init(void)
{
	mm_tbuf_t *b;
	b = (mm_tbuf_t*)calloc(1, sizeof(mm_tbuf_t));
	if (!(mm_dbg_flag & 1)) b->km = km_init();
	return b;
}

void mm_tbuf_destroy(mm_tbuf_t *b)
{
	if (b == 0) return;
	km_destroy(b->km);
	free(b);
}

void *mm_tbuf_get_km(mm_tbuf_t *b)
{
	return b->km;
}

static int mm_dust_minier(void *km, int n, mm128_t *a, int l_seq, const char *seq, int sdust_thres)
{
	int n_dreg, j, k, u = 0;
	const uint64_t *dreg;
	sdust_buf_t *sdb;
	if (sdust_thres <= 0) return n;
	sdb = sdust_buf_init(km);
	dreg = sdust_core((const uint8_t*)seq, l_seq, sdust_thres, 64, &n_dreg, sdb);
	for (j = k = 0; j < n; ++j) { // squeeze out minimizers that significantly overlap with LCRs
		int32_t qpos = (uint32_t)a[j].y>>1, span = a[j].x&0xff;
		int32_t s = qpos - (span - 1), e = s + span;
		while (u < n_dreg && (int32_t)dreg[u] <= s) ++u;
		if (u < n_dreg && (int32_t)(dreg[u]>>32) < e) {
			int v, l = 0;
			for (v = u; v < n_dreg && (int32_t)(dreg[v]>>32) < e; ++v) { // iterate over LCRs overlapping this minimizer
				int ss = s > (int32_t)(dreg[v]>>32)? s : dreg[v]>>32;
				int ee = e < (int32_t)dreg[v]? e : (uint32_t)dreg[v];
				l += ee - ss;
			}
			if (l <= span>>1) a[k++] = a[j]; // keep the minimizer if less than half of it falls in masked region
		} else a[k++] = a[j];
	}
	sdust_buf_destroy(sdb);
	return k; // the new size
}

static void collect_minimizers(void *km, const mm_mapopt_t *opt, const mm_idx_t *mi, int n_segs, const int *qlens, const char **seqs, mm128_v *mv)
{
	int i, n, sum = 0;
	mv->n = 0;
	for (i = n = 0; i < n_segs; ++i) {
		size_t j;
		mm_sketch(km, seqs[i], qlens[i], mi->w, mi->k, i, mi->flag&MM_I_HPC, mv);
		for (j = n; j < mv->n; ++j)
			mv->a[j].y += sum << 1;
		if (opt->sdust_thres > 0) // mask low-complexity minimizers
			mv->n = n + mm_dust_minier(km, mv->n - n, mv->a + n, qlens[i], seqs[i], opt->sdust_thres);
		sum += qlens[i], n = mv->n;
	}
}

#include "ksort.h"
#define heap_lt(a, b) ((a).x > (b).x)
KSORT_INIT(heap, mm128_t, heap_lt)

static inline int skip_seed(int flag, uint64_t r, const mm_seed_t *q, const char *qname, int qlen, const mm_idx_t *mi, int *is_self)
{
	*is_self = 0;
	if (qname && (flag & (MM_F_NO_DIAG|MM_F_NO_DUAL))) {
		const mm_idx_seq_t *s = &mi->seq[r>>32];
		int cmp;
		cmp = strcmp(qname, s->name);
		if ((flag&MM_F_NO_DIAG) && cmp == 0 && (int)s->len == qlen) {
			if ((uint32_t)r>>1 == (q->q_pos>>1)) return 1; // avoid the diagnonal anchors
			if ((r&1) == (q->q_pos&1)) *is_self = 1; // this flag is used to avoid spurious extension on self chain
		}
		if ((flag&MM_F_NO_DUAL) && cmp > 0) // all-vs-all mode: map once
			return 1;
	}
	if (flag & (MM_F_FOR_ONLY|MM_F_REV_ONLY)) {
		if ((r&1) == (q->q_pos&1)) { // forward strand
			if (flag & MM_F_REV_ONLY) return 1;
		} else {
			if (flag & MM_F_FOR_ONLY) return 1;
		}
	}
	return 0;
}

static mm128_t *collect_seed_hits_heap(void *km, const mm_mapopt_t *opt, int max_occ, const mm_idx_t *mi, const char *qname, const mm128_v *mv, int qlen, int64_t *n_a, int *rep_len,
								  int *n_mini_pos, uint64_t **mini_pos)
{
	int i, n_m, heap_size = 0;
	int64_t j, n_for = 0, n_rev = 0;
	mm_seed_t *m;
	mm128_t *a, *heap;

	m = mm_collect_matches(km, &n_m, qlen, max_occ, opt->max_max_occ, opt->occ_dist, mi, mv, n_a, rep_len, n_mini_pos, mini_pos);

	heap = (mm128_t*)kmalloc(km, n_m * sizeof(mm128_t));
	a = (mm128_t*)kmalloc(km, *n_a * sizeof(mm128_t));

	for (i = 0, heap_size = 0; i < n_m; ++i) {
		if (m[i].n > 0) {
			heap[heap_size].x = m[i].cr[0];
			heap[heap_size].y = (uint64_t)i<<32;
			++heap_size;
		}
	}
	ks_heapmake_heap(heap_size, heap);
	while (heap_size > 0) {
		mm_seed_t *q = &m[heap->y>>32];
		mm128_t *p;
		uint64_t r = heap->x;
		int32_t is_self, rpos = (uint32_t)r >> 1;
		if (!skip_seed(opt->flag, r, q, qname, qlen, mi, &is_self)) {
			if ((r&1) == (q->q_pos&1)) { // forward strand
				p = &a[n_for++];
				p->x = (r&0xffffffff00000000ULL) | rpos;
				p->y = (uint64_t)q->q_span << 32 | q->q_pos >> 1;
			} else { // reverse strand
				p = &a[(*n_a) - (++n_rev)];
				p->x = 1ULL<<63 | (r&0xffffffff00000000ULL) | rpos;
				p->y = (uint64_t)q->q_span << 32 | (qlen - ((q->q_pos>>1) + 1 - q->q_span) - 1);
			}
			p->y |= (uint64_t)q->seg_id << MM_SEED_SEG_SHIFT;
			if (q->is_tandem) p->y |= MM_SEED_TANDEM;
			if (is_self) p->y |= MM_SEED_SELF;
		}
		// update the heap
		if ((uint32_t)heap->y < q->n - 1) {
			++heap[0].y;
			heap[0].x = m[heap[0].y>>32].cr[(uint32_t)heap[0].y];
		} else {
			heap[0] = heap[heap_size - 1];
			--heap_size;
		}
		ks_heapdown_heap(0, heap_size, heap);
	}
	kfree(km, m);
	kfree(km, heap);

	// reverse anchors on the reverse strand, as they are in the descending order
	for (j = 0; j < n_rev>>1; ++j) {
		mm128_t t = a[(*n_a) - 1 - j];
		a[(*n_a) - 1 - j] = a[(*n_a) - (n_rev - j)];
		a[(*n_a) - (n_rev - j)] = t;
	}
	if (*n_a > n_for + n_rev) {
		memmove(a + n_for, a + (*n_a) - n_rev, n_rev * sizeof(mm128_t));
		*n_a = n_for + n_rev;
	}
	return a;
}

static mm128_t *collect_seed_hits(void *km, const mm_mapopt_t *opt, int max_occ, const mm_idx_t *mi, const char *qname, const mm128_v *mv, int qlen, int64_t *n_a, int *rep_len,
								  int *n_mini_pos, uint64_t **mini_pos)
{
	int i, n_m;
	mm_seed_t *m;
	mm128_t *a;
	m = mm_collect_matches(km, &n_m, qlen, max_occ, opt->max_max_occ, opt->occ_dist, mi, mv, n_a, rep_len, n_mini_pos, mini_pos);
	a = (mm128_t*)kmalloc(km, *n_a * sizeof(mm128_t));
	for (i = 0, *n_a = 0; i < n_m; ++i) {
		mm_seed_t *q = &m[i];
		const uint64_t *r = q->cr;
		uint32_t k;
		for (k = 0; k < q->n; ++k) {
			int32_t is_self, rpos = (uint32_t)r[k] >> 1;
			mm128_t *p;
			if (skip_seed(opt->flag, r[k], q, qname, qlen, mi, &is_self)) continue;
			p = &a[(*n_a)++];
			if ((r[k]&1) == (q->q_pos&1)) { // forward strand
				p->x = (r[k]&0xffffffff00000000ULL) | rpos;
				p->y = (uint64_t)q->q_span << 32 | q->q_pos >> 1;
			} else if (!(opt->flag & MM_F_QSTRAND)) { // reverse strand and not in the query-strand mode
				p->x = 1ULL<<63 | (r[k]&0xffffffff00000000ULL) | rpos;
				p->y = (uint64_t)q->q_span << 32 | (qlen - ((q->q_pos>>1) + 1 - q->q_span) - 1);
			} else { // reverse strand; query-strand
				int32_t len = mi->seq[r[k]>>32].len;
				p->x = 1ULL<<63 | (r[k]&0xffffffff00000000ULL) | (len - (rpos + 1 - q->q_span) - 1); // coordinate only accurate for non-HPC seeds
				p->y = (uint64_t)q->q_span << 32 | q->q_pos >> 1;
			}
			p->y |= (uint64_t)q->seg_id << MM_SEED_SEG_SHIFT;
			if (q->is_tandem) p->y |= MM_SEED_TANDEM;
			if (is_self) p->y |= MM_SEED_SELF;
		}
	}
	kfree(km, m);
	radix_sort_128x(a, a + (*n_a));
	return a;
}

static void chain_post(const mm_mapopt_t *opt, int max_chain_gap_ref, const mm_idx_t *mi, void *km, int qlen, int n_segs, const int *qlens, int *n_regs, mm_reg1_t *regs, mm128_t *a)
{
	if (!(opt->flag & MM_F_ALL_CHAINS)) { // don't choose primary mapping(s)
		mm_set_parent(km, opt->mask_level, opt->mask_len, *n_regs, regs, opt->a * 2 + opt->b, opt->flag&MM_F_HARD_MLEVEL, opt->alt_drop);
		if (n_segs <= 1) mm_select_sub(km, opt->pri_ratio, mi->k*2, opt->best_n, 1, opt->max_gap * 0.8, n_regs, regs);
		else mm_select_sub_multi(km, opt->pri_ratio, 0.2f, 0.7f, max_chain_gap_ref, mi->k*2, opt->best_n, n_segs, qlens, n_regs, regs);
	}
}

static mm_reg1_t *align_regs(const mm_mapopt_t *opt, const mm_idx_t *mi, void *km, int qlen, const char *seq, int *n_regs, mm_reg1_t *regs, mm128_t *a)
{
	if (!(opt->flag & MM_F_CIGAR)) return regs;
	regs = mm_align_skeleton(km, opt, mi, qlen, seq, n_regs, regs, a); // this calls mm_filter_regs()
	if (!(opt->flag & MM_F_ALL_CHAINS)) { // don't choose primary mapping(s)
		mm_set_parent(km, opt->mask_level, opt->mask_len, *n_regs, regs, opt->a * 2 + opt->b, opt->flag&MM_F_HARD_MLEVEL, opt->alt_drop);
		mm_select_sub(km, opt->pri_ratio, mi->k*2, opt->best_n, 0, opt->max_gap * 0.8, n_regs, regs);
		mm_set_sam_pri(*n_regs, regs);
	}
	return regs;
}

// collecting unique chromosome names found for the current read
/* ------------------------------------------------------------------------
 * Equivalence of candidate alignments by projection onto the base reference.
 *
 * An augmented reference is T' = T + contigs, and the contigs are named
 * "<chrom>_contig_<pos>", where <pos> is the coordinate of the variant on
 * <chrom>. A contig is therefore a copy of T[pos-w .. pos+w] with one variant
 * applied, so *any* alignment inside it projects, under the coordinate map
 * lambda, into [pos - L, pos + L] on <chrom>, where L is the contig length.
 * That bound is valid for every flank w <= L and needs no extra input file:
 * the name and the length are both already in the index.
 *
 * Two candidates are equivalent iff their projections overlap as intervals
 * in T. This is the relation ~ from the formal model; the name-only rule it
 * replaces discards every contig candidate regardless of locus, which merges
 * candidates that really are independent hypotheses and inflates mapq.
 *
 * Switched with MMQ_PROJECTION=1 so that one binary can run both rules.
 * Before publication this belongs in options.c as a proper flag.
 * ---------------------------------------------------------------------- */
typedef struct {
	const char *chrom;
	size_t chrom_len;
	int64_t st, en;
} mmq_proj_t;

static void mmq_project(const mm_idx_t *mi, const mm_reg1_t *r, mmq_proj_t *pr)
{
	const char *name = mi->seq[r->rid].name;
	const char *tag  = strstr(name, "_contig_");
	if (tag == NULL) {                       /* already a sequence of T */
		pr->chrom = name;
		pr->chrom_len = strlen(name);
		pr->st = r->rs;
		pr->en = r->re;
	} else {
		int64_t pos = strtoll(tag + 8, NULL, 10);
		int64_t len = (int64_t)mi->seq[r->rid].len;
		/* контиг = T[p-w .. p-1] + alt + T[p+|ref| .. p+w-1], значит |c| = 2w + |alt| >= 2w,
		 * откуда w <= |c|/2, а проекция лежит в [p - |c|/2, p + |c|/2].
		 * Это по-прежнему оценка сверху, но вдвое точнее прежней [p-|c|, p+|c|]. */
		int64_t half = len / 2 + 1;
		pr->chrom = name;
		pr->chrom_len = (size_t)(tag - name);
		pr->st = pos - half;
		pr->en = pos + half;
	}
}

static int mmq_equivalent(const mmq_proj_t *a, const mmq_proj_t *b)
{
	if (a->chrom_len != b->chrom_len) return 0;
	if (strncmp(a->chrom, b->chrom, a->chrom_len) != 0) return 0;
	return a->st <= b->en && b->st <= a->en;
}

/* Same-locus test on two candidates.  With --lift the candidates have already
 * been projected exactly (mm_lift_regs runs BEFORE this rule), so the test is
 * mm_lift_same_locus(): same sequence, same strand, a shared diagonal.  Without
 * the projection table the old name-parsed upper bound is all there is. */
static int mmq_same_locus(const mm_idx_t *mi, const mm_reg1_t *a, const mm_reg1_t *b)
{
	mmq_proj_t pa, pb;
	if (mi->lift) return mm_lift_same_locus(a, b);
	mmq_project(mi, a, &pa);
	mmq_project(mi, b, &pb);
	return mmq_equivalent(&pa, &pb);
}

/* "Came from a contig": the name before projection, the lifted bit after. */
static int mmq_is_contig(const mm_idx_t *mi, const mm_reg1_t *r)
{
	if (mi->lift) return r->lifted;
	return strstr(mi->seq[r->rid].name, "_contig_") != NULL;
}


/* Disagreement between the two rules, counted on whatever data is being mapped.
 * mmq_name_only is the case the paper cares about: the name rule discards a
 * candidate whose projection does NOT meet the best one's, i.e. it destroys a
 * genuinely independent hypothesis. Reported at exit; costs one comparison
 * per candidate, so it can be left on for a production run. */
#define MMQ_RULE(v) __atomic_fetch_add(&(v), 1, __ATOMIC_RELAXED)
static long mmq_agree = 0, mmq_name_only = 0, mmq_proj_only = 0, mmq_neither = 0;
__attribute__((destructor)) static void mmq_rule_report(void)
{
	long tot = mmq_agree + mmq_name_only + mmq_proj_only + mmq_neither;
	if (tot == 0) return;
	fprintf(stderr, "[mmq-rules] candidates=%ld agree=%ld name_only=%ld proj_only=%ld neither=%ld"
			" | name_only/all=%.3f%%",
			tot, mmq_agree, mmq_name_only, mmq_proj_only, mmq_neither,
			tot ? 100.0 * mmq_name_only / tot : 0.0);
	{ long coll = mmq_agree + mmq_name_only;
	  fprintf(stderr, " name_only/collapsed=%.1f%%\n", coll ? 100.0 * mmq_name_only / coll : 0.0); }
}


/* Among equally-scoring members of the best candidate's class, keep the one
 * that already lives on T. lambda is the identity there, so the read needs no
 * lift-over, its mate stays in the same sequence, and its CIGAR is not split
 * by the coordinate map. Only exact score ties are moved, so the alignment
 * decision itself is untouched. */
static void mmq_prefer_t_representative(const mm_idx_t *mi, int n_a, mm_reg1_t *regs)
{
	mmq_proj_t p0, pi;
	int i;
	(void)p0; (void)pi;
	if (n_a < 2 || regs == NULL) return;
	if (!mmq_is_contig(mi, &regs[0])) return; /* already on T */
	if (regs[0].p == NULL) return;
	for (i = 1; i < n_a; ++i) {
		if (regs[i].p == NULL) continue;
		if (mmq_is_contig(mi, &regs[i])) continue;
		if (regs[i].p->dp_score != regs[0].p->dp_score) continue;
		if (!mmq_same_locus(mi, &regs[0], &regs[i])) continue;
		{
			/* sam_pri, id and parent belong to the slot, not to the record:
			 * slot 0 stays the primary, only the alignment itself moves. */
			mm_reg1_t tmp = regs[0];
			uint32_t pri0 = regs[0].sam_pri, prii = regs[i].sam_pri;
			int id0 = regs[0].id, idi = regs[i].id;
			int par0 = regs[0].parent, pari = regs[i].parent;
			regs[0] = regs[i];      regs[i] = tmp;
			regs[0].id = id0;       regs[i].id = idi;
			regs[0].parent = par0;  regs[i].parent = pari;
			regs[0].sam_pri = pri0; regs[i].sam_pri = prii;
		}
		return;
	}
}

char **collect_seed_chromosome_names(const mm_idx_t *mi, int n_a, mm_reg1_t *regs, int64_t flag)
{
	char **chromosome_names = CALLOC(char*, n_a);
	const char* best_alignment;
	char* best_copy;
	const char* delete_name = "delete";
	int i;
	for (i = 0; i < n_a; ++i) {
		chromosome_names[i] = NULL;
	}
	if (regs != NULL) {
		if (regs[0].parent != 0) {
			int parent_id = regs[0].parent;
			mm_reg1_t tmp = regs[0];
			regs[0] = regs[parent_id];
			regs[parent_id] = tmp;

			regs[parent_id].id = regs[0].id;
			regs[0].id = tmp.id;

			regs[0].parent = 0;
			regs[parent_id].parent = 0;
		}

		/* сделать до цикла: перестановка меняет regs[0], а имена читаются внутри цикла */
		if (flag & MM_F_MAPQ_PREFER_T) mmq_prefer_t_representative(mi, n_a, regs);

		for (i = 1; i < n_a; ++i) {
			const char *curr_name;
			char *curr_copy;
			char *output;
			char *best;
			regs[i].parent = 0;
			best_alignment = mi->seq[regs[0].rid].name;
			best_copy = strdup(best_alignment);
			curr_name = mi->seq[regs[i].rid].name;
			curr_copy = strdup(curr_name);

			{	/* count how often the two rules disagree, whichever is in force */
				int by_name, by_proj;
				by_proj = mmq_same_locus(mi, &regs[0], &regs[i]);
				by_name = mmq_is_contig(mi, &regs[i]);
				if (!by_name && mmq_is_contig(mi, &regs[0])) {
					char *b2 = strdup(best_alignment), *c2 = strdup(curr_name);
					char *bt = strtok(b2, "_"), *ct = strtok(c2, "_");
					by_name = (bt && ct && strcmp(bt, ct) == 0);
					free(b2); free(c2);
				}
				if (by_name && by_proj)        MMQ_RULE(mmq_agree);
				else if (by_name && !by_proj)  MMQ_RULE(mmq_name_only);
				else if (!by_name && by_proj)  MMQ_RULE(mmq_proj_only);
				else                           MMQ_RULE(mmq_neither);
			}
			if (flag & MM_F_MAPQ_QUOTIENT) {
				/* collapse only candidates that are the same locus as the best:
				 * symmetric -- it does not matter which of the two is the contig */
				if (mmq_same_locus(mi, &regs[0], &regs[i]))
					chromosome_names[i] = strdup(delete_name);
			} else {
			if (strstr(curr_name, "contig") != NULL) {
				chromosome_names[i] = strdup(delete_name);
			}
			if (strstr(best_alignment, "contig") != NULL) {
				best = strtok(best_copy, "_");
				output = strtok(curr_copy, "_");
				if (best != NULL && output != NULL && strcmp(best, output) == 0) {
					chromosome_names[i] = strdup(delete_name);
				}
			}
			}

			free(curr_copy);
			free(best_copy);
		}
	}

	return chromosome_names;
}

void free_chromosome_names(char **chromosome_names, int max_chromosomes) {
	int i;
	for (i = 0; i < max_chromosomes; i++) {
		if (chromosome_names[i] != NULL) {
			free(chromosome_names[i]);
		}
	}
	free(chromosome_names);
}

typedef struct {
	char *name;			// chromosome name
	uint32_t len;		// chromosome length
	int32_t rid;		// chromosome id
	int32_t n_hits;		// number of hits to the current chromosome
	int32_t start_pos;	// start position of the first hit
	int32_t end_pos;	// end position of the last hit
} chromosome_info_t;

// collecting unique chromosome info found for the current read
chromosome_info_t *collect_seed_chromosome_info(const mm_idx_t *mi, int n_a, mm128_t *a, int *n_chromosomes) 
{
	int max_chroms = 32;
	*n_chromosomes = 0;
	chromosome_info_t *chr_info = (chromosome_info_t*)calloc(max_chroms, sizeof(chromosome_info_t));
	int i;

	for (i = 0; i < n_a; ++i) {
		uint32_t rid = a[i].x<<1>>33;
		int exists = 0;
		int32_t pos = (int32_t)a[i].x;
		int j;

		for (j = 0; j < *n_chromosomes; ++j) {
			if (rid == chr_info[j].rid) {
				exists = 1;
				chr_info[j].n_hits++;

				if (pos < chr_info[j].start_pos) chr_info[j].start_pos = pos;
				if (pos > chr_info[j].end_pos) chr_info[j].end_pos = pos;

				break;
			}
		}

		if (!exists) {
			if (*n_chromosomes == max_chroms) {
				max_chroms *= 2;
				chr_info = (chromosome_info_t*)realloc(chr_info, max_chroms * sizeof(chromosome_info_t));
			}

			int idx = *n_chromosomes;
			chr_info[idx].name = strdup(mi->seq[rid].name);
			chr_info[idx].len = mi->seq[rid].len;
			chr_info[idx].rid = rid;
			chr_info[idx].n_hits = 1;
			chr_info[idx].start_pos = pos;
			chr_info[idx].end_pos = pos;
			(*n_chromosomes)++;
		}
	}

	if (*n_chromosomes < max_chroms) {
		chr_info = (chromosome_info_t*)realloc(chr_info, *n_chromosomes * sizeof(chromosome_info_t));
	}

	return chr_info;
}

// free chromosome_info
void free_chromosome_info(chromosome_info_t *chr_info, int n_chromosomes) {
	int i;
	for (i = 0; i < n_chromosomes; i++) {
		free(chr_info[i].name);
	}
	free(chr_info);
}

// delete second alignment (to improve mapq)
mm_reg1_t* remove_second_suboptimal_alignment(mm_reg1_t *regs, int *num_regs, int z, int sub_diff) {
	if (!regs || *num_regs < 2) {
		return regs;
	}

	if (regs[z].p) {
		free(regs[z].p);
	}
	int i;
	for (i = z; i < *num_regs - 1; i++) {
		regs[i] = regs[i + 1];
	}

	(*num_regs)--;

	regs = (mm_reg1_t*)realloc(regs, (*num_regs) * sizeof(*regs));

	return regs;
}

/* The fork's per-segment sequence: locus rule -> drop collapsed candidates ->
 * refresh dp_max2/subsc (and n_sub when asked) -> mapq.  Factored out so the
 * mate step can re-run it on a pruned candidate set. */
static mm_reg1_t *mmq_rule_and_mapq(const mm_idx_t *mi, const mm_mapopt_t *opt, void *km, int *n_regs_, mm_reg1_t *regs,
                                    int rep_len, int is_sr, int is_sr_rna, int is_splice, int recompute_nsub, int apply_rule)
{
	int n_regs = *n_regs_, z, n_names = n_regs;
	char **chrs_to_drop = NULL;
	const char *delete_name = "delete";
	if (apply_rule) {
		chrs_to_drop = collect_seed_chromosome_names(mi, n_regs, regs, opt->flag);
		for (z = n_regs - 1; z > 0; z--)
			if (chrs_to_drop[z] != NULL && strcmp(chrs_to_drop[z], delete_name) == 0)
				regs = remove_second_suboptimal_alignment(regs, &n_regs, z, opt->a * 2 + opt->b);
	}
	if (n_regs > 0) {
		int max_dpmax2 = 0, max_score = 0;
		for (z = 1; z < n_regs; z++) {
			regs[z].id = z;
			if (regs[z].p->dp_score > max_dpmax2) max_dpmax2 = regs[z].p->dp_score;
			if (regs[z].score > max_score) max_score = regs[z].score;
		}
		if (regs[0].p != NULL) regs[0].p->dp_max2 = max_dpmax2;
		regs[0].subsc = max_score;
		if (recompute_nsub) {
			mm_reg1_t *copy = (mm_reg1_t*)malloc(n_regs * sizeof(mm_reg1_t));
			memcpy(copy, regs, n_regs * sizeof(mm_reg1_t));
			copy[0].n_sub = 0;
			mm_set_parent(km, opt->mask_level, opt->mask_len, n_regs, copy, opt->a * 2 + opt->b, opt->flag&MM_F_HARD_MLEVEL, opt->alt_drop);
			regs[0].n_sub = copy[0].n_sub;
			free(copy);
		}
	}
	mm_set_mapq2(km, n_regs, regs, opt->min_chain_score, opt->a, rep_len, is_sr || is_sr_rna, is_splice, chrs_to_drop);
	if (chrs_to_drop) free_chromosome_names(chrs_to_drop, n_names);
	*n_regs_ = n_regs;
	return regs;
}

/* A candidate c is compatible with the (unambiguous) mate m when both lie on
 * one sequence, on opposite strands (FR), and the fragment they imply is no
 * longer than max_frag_len.  Returns the implied fragment length, or -1. */
static int mm_mate_compat(const mm_mapopt_t *opt, int s, const mm_reg1_t *c, const mm_reg1_t *m)
{
	int lo, hi, span;
	/* In fragment mode the second read is reverse-complemented before mapping
	 * (mm_map_frag), so a proper FR pair has BOTH segments with the same rev;
	 * mm_pair() encodes the same convention (key bit = s ^ rev). */
	if (c->p == NULL || c->rid != m->rid || c->rev != m->rev) return -1;
	if (opt->pe_ori == 1) {	/* FR: the left end is segment 0 when rev == 0, segment 1 when rev == 1 */
		const mm_reg1_t *seg0 = s == 0 ? c : m, *seg1 = s == 0 ? m : c;
		const mm_reg1_t *left = c->rev ? seg1 : seg0, *right = c->rev ? seg0 : seg1;
		if (left->rs > right->re) return -1;
	}
	lo = c->rs < m->rs ? c->rs : m->rs;
	hi = c->re > m->re ? c->re : m->re;
	span = hi - lo;
	if (span > opt->max_frag_len) return -1;
	return span;
}

/* --mate-resolve, BEFORE mapq is final and BEFORE mm_pair(): if exactly one
 * end is unambiguous (standalone mapq >= opt->mate_mapq), the other end keeps
 * only the candidate nearest to that mate and its mapq is recomputed on that
 * set.  If both ends are ambiguous nothing is touched and the pair is marked
 * (pa:i:1) so those cases can be pulled out and looked at.  Returns 1 when
 * both ends were ambiguous. */
static int mm_mate_resolve(const mm_idx_t *mi, const mm_mapopt_t *opt, void *km, int *n_regs, mm_reg1_t **regs,
                           int rep_len, int is_sr, int is_sr_rna, int is_splice)
{
	int s, t, i, best = -1, best_span = INT32_MAX, n_compat = 0, unamb[2], maxdp = 0;
	mm_reg1_t *anchor;
	if (n_regs[0] == 0 || n_regs[1] == 0) return 0;
	/* "unambiguous" = standalone mapq at least --mate-mapq AND no competitor
	 * within sub_diff of the best DP score (n_sub == 0).  The second condition
	 * matters because for short paired reads the fragment chain that includes
	 * the mate inflates the primary's chaining score, and the mapq formula then
	 * reports 60 for a read with six equally-scoring copies. */
	for (s = 0; s < 2; ++s) unamb[s] = regs[s][0].mapq >= opt->mate_mapq && regs[s][0].n_sub == 0;
	if (mm_dbg_flag & MM_DBG_PRINT_QNAME)
		fprintf(stderr, "MR0\tn=%d,%d slot0 mapq=%d,%d id/parent=%d/%d,%d/%d sam_pri=%d,%d\n", n_regs[0], n_regs[1], regs[0][0].mapq, regs[1][0].mapq,
				regs[0][0].id, regs[0][0].parent, regs[1][0].id, regs[1][0].parent, regs[0][0].sam_pri, regs[1][0].sam_pri);
	if (mm_dbg_flag & MM_DBG_PRINT_QNAME)
		for (s = 0; s < 2; ++s) for (i = 0; i < n_regs[s]; ++i)
			fprintf(stderr, "MR0\t  seg%d[%d] %d:%d-%d q%d-%d rev=%d id=%d par=%d score=%d dp=%d dp_max2=%d subsc=%d n_sub=%d mapq=%d lifted=%d\n", s, i,
					regs[s][i].rid, regs[s][i].rs, regs[s][i].re, regs[s][i].qs, regs[s][i].qe, regs[s][i].rev, regs[s][i].id, regs[s][i].parent, regs[s][i].score,
					regs[s][i].p ? regs[s][i].p->dp_score : -1, regs[s][i].p ? regs[s][i].p->dp_max2 : -1, regs[s][i].subsc, regs[s][i].n_sub, regs[s][i].mapq, regs[s][i].lifted);
	if (unamb[0] && unamb[1]) return 0;
	if (!unamb[0] && !unamb[1]) return 1;
	s = unamb[0] ? 1 : 0; t = 1 - s;
	anchor = &regs[t][0];
	/* only the candidates that made the end ambiguous take part: those within
	 * sub_diff (one mismatch) of the best DP score -- the same window
	 * mm_set_parent() uses for n_sub.  A clipped partial alignment that happens
	 * to sit nearest the mate is not a hypothesis the read was torn between. */
	for (i = 0; i < n_regs[s]; ++i)
		if (regs[s][i].p && regs[s][i].p->dp_max > maxdp) maxdp = regs[s][i].p->dp_max;
	if (mm_dbg_flag & MM_DBG_PRINT_QNAME) fprintf(stderr, "MR\tseg=%d n=%d mapq=%d anchor=%d:%d-%d rev=%d\n", s, n_regs[s], regs[s][0].mapq, anchor->rid, anchor->rs, anchor->re, anchor->rev);
	for (i = 0; i < n_regs[s]; ++i) {
		int span;
		if (regs[s][i].p == NULL || regs[s][i].p->dp_max < maxdp - (opt->a * 2 + opt->b)) continue;
		span = mm_mate_compat(opt, s, &regs[s][i], anchor);
		if (mm_dbg_flag & MM_DBG_PRINT_QNAME) fprintf(stderr, "MR\t  cand %d: %d:%d-%d q%d-%d rev=%d dp_max=%d dp_score=%d span=%d\n", i, regs[s][i].rid, regs[s][i].rs, regs[s][i].re, regs[s][i].qs, regs[s][i].qe, regs[s][i].rev, regs[s][i].p->dp_max, regs[s][i].p->dp_score, span);
		if (span < 0) continue;
		++n_compat;
		if (span < best_span) best_span = span, best = i;   /* ties: lower index = better score */
	}
	if (best < 0) return 0;   /* nothing compatible: leave it to mm_pair() as before */
	/* drop everything but the nearest compatible candidate, highest index first */
	for (i = n_regs[s] - 1; i >= 0 && n_regs[s] > 1; --i) {
		if (i == best) continue;
		regs[s] = remove_second_suboptimal_alignment(regs[s], &n_regs[s], i, opt->a * 2 + opt->b);
		if (i < best) --best;
	}
	/* slot 0 is the primary; the survivor is there now (best == 0) */
	regs[s][0].id = 0; regs[s][0].parent = 0; regs[s][0].sam_pri = 1;
	regs[s][0].mate_ncand = n_compat; /* emitted as mr:i:N when N > 1 -- the "nearest of several" cases to look at by eye */
	regs[s] = mmq_rule_and_mapq(mi, opt, km, &n_regs[s], regs[s], rep_len, is_sr, is_sr_rna, is_splice, 1, 0);
	return 0;
}

void mm_map_frag_core(const mm_idx_t *mi, int n_segs, const int *qlens, const char **seqs, int *n_regs, mm_reg1_t **regs, mm_tbuf_t *b, const mm_mapopt_t *opt, const char *qname)
{
	int i, j, rep_len, qlen_sum, n_regs0, n_mini_pos;
	int max_chain_gap_qry, max_chain_gap_ref, is_splice = !!(opt->flag & MM_F_SPLICE), is_sr = !!(opt->flag & MM_F_SR), is_sr_rna = !!(opt->flag & MM_F_SR_RNA);
	uint32_t hash;
	int64_t n_a;
	uint64_t *u, *mini_pos;
	mm128_t *a;
	mm128_v mv = {0,0,0};
	mm_reg1_t *regs0;
	km_stat_t kmst;
	float chn_pen_gap, chn_pen_skip;

	for (i = 0, qlen_sum = 0; i < n_segs; ++i)
		qlen_sum += qlens[i], n_regs[i] = 0, regs[i] = 0;

	if (qlen_sum == 0 || n_segs <= 0 || n_segs > MM_MAX_SEG) return;
	if (opt->max_qlen > 0 && qlen_sum > opt->max_qlen) return;

	hash  = qname && !(opt->flag & MM_F_NO_HASH_NAME)? __ac_X31_hash_string(qname) : 0;
	hash ^= __ac_Wang_hash(qlen_sum) + __ac_Wang_hash(opt->seed);
	hash  = __ac_Wang_hash(hash);

	collect_minimizers(b->km, opt, mi, n_segs, qlens, seqs, &mv);
	if (opt->q_occ_frac > 0.0f) mm_seed_mz_flt(b->km, &mv, opt->mid_occ, opt->q_occ_frac);
	if (opt->flag & MM_F_HEAP_SORT) a = collect_seed_hits_heap(b->km, opt, opt->mid_occ, mi, qname, &mv, qlen_sum, &n_a, &rep_len, &n_mini_pos, &mini_pos);
	else a = collect_seed_hits(b->km, opt, opt->mid_occ, mi, qname, &mv, qlen_sum, &n_a, &rep_len, &n_mini_pos, &mini_pos);

	if (mm_dbg_flag & MM_DBG_PRINT_SEED) {
		fprintf(stderr, "RS\t%d\n", rep_len);
		for (i = 0; i < n_a; ++i)
			fprintf(stderr, "SD\t%s\t%d\t%c\t%d\t%d\t%d\n", mi->seq[a[i].x<<1>>33].name, (int32_t)a[i].x, "+-"[a[i].x>>63], (int32_t)a[i].y, (int32_t)(a[i].y>>32&0xff),
					i == 0? 0 : ((int32_t)a[i].y - (int32_t)a[i-1].y) - ((int32_t)a[i].x - (int32_t)a[i-1].x));
	}

	// set max chaining gap on the query and the reference sequence
	if (is_sr)
		max_chain_gap_qry = qlen_sum > opt->max_gap? qlen_sum : opt->max_gap;
	else max_chain_gap_qry = opt->max_gap;
	if (opt->max_gap_ref > 0) {
		max_chain_gap_ref = opt->max_gap_ref; // always honor mm_mapopt_t::max_gap_ref if set
	} else if (opt->max_frag_len > 0) {
		max_chain_gap_ref = opt->max_frag_len - qlen_sum;
		if (max_chain_gap_ref < opt->max_gap) max_chain_gap_ref = opt->max_gap;
	} else max_chain_gap_ref = opt->max_gap;

	chn_pen_gap  = opt->chain_gap_scale * 0.01 * mi->k;
	chn_pen_skip = opt->chain_skip_scale * 0.01 * mi->k;
	if (opt->flag & MM_F_RMQ) {
		a = mg_lchain_rmq(opt->max_gap, opt->rmq_inner_dist, opt->bw, opt->max_chain_skip, opt->rmq_size_cap, opt->min_cnt, opt->min_chain_score,
						  chn_pen_gap, chn_pen_skip, n_a, a, &n_regs0, &u, b->km);
	} else {
		a = mg_lchain_dp(max_chain_gap_ref, max_chain_gap_qry, opt->bw, opt->max_chain_skip, opt->max_chain_iter, opt->min_cnt, opt->min_chain_score,
						 chn_pen_gap, chn_pen_skip, is_splice, n_segs, n_a, a, &n_regs0, &u, b->km);
	}

	if (opt->bw_long > opt->bw && (opt->flag & (MM_F_SPLICE|MM_F_SR|MM_F_NO_LJOIN)) == 0 && n_segs == 1 && n_regs0 > 1) { // re-chain/long-join for long sequences
		int32_t st = (int32_t)a[0].y, en = (int32_t)a[(int32_t)u[0] - 1].y;
		if (qlen_sum - (en - st) > opt->rmq_rescue_size || en - st > qlen_sum * opt->rmq_rescue_ratio) {
			int32_t i;
			for (i = 0, n_a = 0; i < n_regs0; ++i) n_a += (int32_t)u[i];
			kfree(b->km, u);
			radix_sort_128x(a, a + n_a);
			a = mg_lchain_rmq(opt->max_gap, opt->rmq_inner_dist, opt->bw_long, opt->max_chain_skip, opt->rmq_size_cap, opt->min_cnt, opt->min_chain_score,
							  chn_pen_gap, chn_pen_skip, n_a, a, &n_regs0, &u, b->km);
		}
	} else if (opt->max_occ > opt->mid_occ && rep_len > 0 && !(opt->flag & MM_F_RMQ)) { // re-chain, mostly for short reads
		int rechain = 0;
		if (n_regs0 > 0) { // test if the best chain has all the segments
			int n_chained_segs = 1, max = 0, max_i = -1, max_off = -1, off = 0;
			for (i = 0; i < n_regs0; ++i) { // find the best chain
				if (max < (int)(u[i]>>32)) max = u[i]>>32, max_i = i, max_off = off;
				off += (uint32_t)u[i];
			}
			for (i = 1; i < (int32_t)u[max_i]; ++i) // count the number of segments in the best chain
				if ((a[max_off+i].y&MM_SEED_SEG_MASK) != (a[max_off+i-1].y&MM_SEED_SEG_MASK))
					++n_chained_segs;
			if (n_chained_segs < n_segs)
				rechain = 1;
		} else rechain = 1;
		if (rechain) { // redo chaining with a higher max_occ threshold
			kfree(b->km, a);
			kfree(b->km, u);
			kfree(b->km, mini_pos);
			if (opt->flag & MM_F_HEAP_SORT) a = collect_seed_hits_heap(b->km, opt, opt->max_occ, mi, qname, &mv, qlen_sum, &n_a, &rep_len, &n_mini_pos, &mini_pos);
			else a = collect_seed_hits(b->km, opt, opt->max_occ, mi, qname, &mv, qlen_sum, &n_a, &rep_len, &n_mini_pos, &mini_pos);
			a = mg_lchain_dp(max_chain_gap_ref, max_chain_gap_qry, opt->bw, opt->max_chain_skip, opt->max_chain_iter, opt->min_cnt, opt->min_chain_score,
							 chn_pen_gap, chn_pen_skip, is_splice, n_segs, n_a, a, &n_regs0, &u, b->km);
		}
	}
	b->frag_gap = max_chain_gap_ref;
	b->rep_len = rep_len;

	regs0 = mm_gen_regs(b->km, hash, qlen_sum, n_regs0, u, a, !!(opt->flag&MM_F_QSTRAND));
	if ((opt->flag & MM_F_LIFT) && mi->lift) { /* chains on augmentation contigs survive fragment-level selection (pe.c) */
		const mm_lift_t *lt = (const mm_lift_t*)mi->lift;
		for (i = 0; i < n_regs0; ++i)
			regs0[i].on_ctg = (regs0[i].rid >= 0 && regs0[i].rid < lt->n && lt->a[regs0[i].rid].base_rid >= 0);
	}
	if (mi->n_alt) {
		mm_mark_alt(mi, n_regs0, regs0);
		mm_hit_sort(b->km, &n_regs0, regs0, opt->alt_drop); // this step can be merged into mm_gen_regs(); will do if this shows up in profile
	}

	if (mm_dbg_flag & (MM_DBG_PRINT_SEED|MM_DBG_PRINT_CHAIN))
		for (j = 0; j < n_regs0; ++j)
			for (i = regs0[j].as; i < regs0[j].as + regs0[j].cnt; ++i)
				fprintf(stderr, "CN\t%d\t%s\t%d\t%c\t%d\t%d\t%d\n", j, mi->seq[a[i].x<<1>>33].name, (int32_t)a[i].x, "+-"[a[i].x>>63], (int32_t)a[i].y, (int32_t)(a[i].y>>32&0xff),
						i == regs0[j].as? 0 : ((int32_t)a[i].y - (int32_t)a[i-1].y) - ((int32_t)a[i].x - (int32_t)a[i-1].x));

	chain_post(opt, max_chain_gap_ref, mi, b->km, qlen_sum, n_segs, qlens, &n_regs0, regs0, a);
	if (!is_sr && !(opt->flag&MM_F_QSTRAND)) {
		mm_est_err(mi, qlen_sum, n_regs0, regs0, a, n_mini_pos, mini_pos);
		n_regs0 = mm_filter_strand_retained(n_regs0, regs0);
	}

	int do_lift = (opt->flag & MM_F_LIFT) && mi->lift;
	int do_shift = do_lift && !(opt->flag & MM_F_LIFT_NOSHIFT);

	if (n_segs == 1) { // uni-segment
		regs0 = align_regs(opt, mi, b->km, qlens[0], seqs[0], &n_regs0, regs0, a);
		regs0 = (mm_reg1_t*)realloc(regs0, sizeof(*regs0) * n_regs0);
		/* project contig placements onto the base reference BEFORE the locus
		 * rule and mapq, so the rule compares exact base coordinates */
		if (do_lift && mm_lift_regs(mi, n_regs0, regs0) > 0 && do_shift)
			for (i = 0; i < n_regs0; ++i) if (regs0[i].lifted) mm_lift_refresh(mi, qlens[0], seqs[0], &regs0[i], 1);
		regs0 = mmq_rule_and_mapq(mi, opt, b->km, &n_regs0, regs0, rep_len, is_sr, is_sr_rna, is_splice, 0, !(opt->flag & MM_F_NO_MAPQ_RULE));
		if (do_shift) /* NM against the base, after mapq (which used the contig identity) */
			for (i = 0; i < n_regs0; ++i) if (regs0[i].lifted) mm_lift_refresh(mi, qlens[0], seqs[0], &regs0[i], 2);
		n_regs[0] = n_regs0, regs[0] = regs0;
	} else { // multi-segment
		mm_seg_t *seg;
		int both_ambig = 0;
		seg = mm_seg_gen(b->km, hash, n_segs, qlens, n_regs0, regs0, n_regs, regs, a); // split fragment chain to separate segment chains
		free (regs0);
		for (i = 0; i < n_segs; ++i) {
			mm_set_parent(b->km, opt->mask_level, opt->mask_len, n_regs[i], regs[i], opt->a * 2 + opt->b, opt->flag&MM_F_HARD_MLEVEL, opt->alt_drop); // update mm_reg1_t::parent
			regs[i] = align_regs(opt, mi, b->km, qlens[i], seqs[i], &n_regs[i], regs[i], seg[i].a);
			/* project BEFORE the locus rule (exact coordinates) and BEFORE
			 * pairing (a read on a contig and its mate on the chromosome carry
			 * different rids until this runs) */
			if (do_lift && mm_lift_regs(mi, n_regs[i], regs[i]) > 0 && do_shift) {
				int z;
				for (z = 0; z < n_regs[i]; ++z) if (regs[i][z].lifted) mm_lift_refresh(mi, qlens[i], seqs[i], &regs[i][z], 1);
			}
		}
		for (i = 0; i < n_segs; ++i)
			regs[i] = mmq_rule_and_mapq(mi, opt, b->km, &n_regs[i], regs[i], rep_len, is_sr, is_sr_rna, is_splice, 1, !(opt->flag & MM_F_NO_MAPQ_RULE));
		mm_seg_free(b->km, n_segs, seg);
		if (n_segs == 2 && opt->pe_ori >= 0 && (opt->flag&MM_F_CIGAR)) {
			if (opt->flag & MM_F_MATE_RESOLVE)
				both_ambig = mm_mate_resolve(mi, opt, b->km, n_regs, regs, rep_len, is_sr, is_sr_rna, is_splice);
			if (do_shift) /* NM against the base, after every mapq computation (they used the contig identity) */
				for (i = 0; i < n_segs; ++i) { int z; for (z = 0; z < n_regs[i]; ++z) if (regs[i][z].lifted) mm_lift_refresh(mi, qlens[i], seqs[i], &regs[i][z], 2); }
			mm_pair(b->km, max_chain_gap_ref, opt->pe_bonus, opt->a * 2 + opt->b, opt->a, qlens, n_regs, regs, !!(opt->flag & MM_F_MATE_RESOLVE)); // pairing
			if (both_ambig)
				for (i = 0; i < 2; ++i) { int z; for (z = 0; z < n_regs[i]; ++z) if (regs[i][z].sam_pri) regs[i][z].pe_ambig = 1; }
		} else if (do_shift) {
			for (i = 0; i < n_segs; ++i) { int z; for (z = 0; z < n_regs[i]; ++z) if (regs[i][z].lifted) mm_lift_refresh(mi, qlens[i], seqs[i], &regs[i][z], 2); }
		}
	}


	kfree(b->km, mv.a);
	kfree(b->km, a);
	kfree(b->km, u);
	kfree(b->km, mini_pos);

	if (mi->J && n_segs == 1 && is_splice)
		for (i = 0; i < n_regs0; ++i)
			mm_jump_split(b->km, mi, opt, qlens[0], (const uint8_t*)seqs[0], &regs0[i], 0);

	if (b->km) {
		km_stat(b->km, &kmst);
		if (mm_dbg_flag & MM_DBG_PRINT_QNAME)
			fprintf(stderr, "QM\t%s\t%d\tcap=%ld,nCore=%ld,largest=%ld\n", qname, qlen_sum, kmst.capacity, kmst.n_cores, kmst.largest);
		assert(kmst.n_blocks == kmst.n_cores); // otherwise, there is a memory leak
		if (kmst.largest > 1U<<28 || (opt->cap_kalloc > 0 && kmst.capacity > opt->cap_kalloc)) {
			if (mm_dbg_flag & MM_DBG_PRINT_QNAME)
				fprintf(stderr, "[W::%s] reset thread-local memory after read %s\n", __func__, qname);
			km_destroy(b->km);
			b->km = km_init();
		}
	}
}

void mm_map_frag(const mm_idx_t *mi, int n_segs, const int *qlens, const char **seqs, int *n_regs, mm_reg1_t **regs, mm_tbuf_t *b, const mm_mapopt_t *opt, const char *qname)
{
	if ((opt->flag & MM_F_WEAK_PAIRING) && n_segs == 2 && opt->pe_ori >= 0 && (opt->flag&MM_F_CIGAR)) {
		int i;
		for (i = 0; i < n_segs; ++i)
			mm_map_frag_core(mi, 1, &qlens[i], &seqs[i], &n_regs[i], &regs[i], b, opt, qname);
		mm_pair(b->km, opt->max_gap_ref, opt->pe_bonus, opt->a * 2 + opt->b, opt->a, qlens, n_regs, regs, !!(opt->flag & MM_F_MATE_RESOLVE));
	} else {
		mm_map_frag_core(mi, n_segs, qlens, seqs, n_regs, regs, b, opt, qname);
	}
}

mm_reg1_t *mm_map(const mm_idx_t *mi, int qlen, const char *seq, int *n_regs, mm_tbuf_t *b, const mm_mapopt_t *opt, const char *qname)
{
	mm_reg1_t *regs;
	mm_map_frag(mi, 1, &qlen, &seq, n_regs, &regs, b, opt, qname);
	return regs;
}

/**************************
 * Multi-threaded mapping *
 **************************/

typedef struct {
	int n_processed, n_threads, n_fp;
	int64_t mini_batch_size;
	const mm_mapopt_t *opt;
	mm_bseq_file_t **fp;
	const mm_idx_t *mi;
	kstring_t str;

	int n_parts;
	uint32_t *rid_shift;
	FILE *fp_split, **fp_parts;
} pipeline_t;

typedef struct {
	const pipeline_t *p;
	int n_seq, n_frag;
	mm_bseq1_t *seq;
	int *n_reg, *seg_off, *n_seg, *rep_len, *frag_gap;
	mm_reg1_t **reg;
	mm_tbuf_t **buf;
} step_t;

static void worker_for(void *_data, long i, int tid) // kt_for() callback
{
	step_t *s = (step_t*)_data;
	int qlens[MM_MAX_SEG], j, off = s->seg_off[i], pe_ori = s->p->opt->pe_ori;
	const char *qseqs[MM_MAX_SEG];
	double t = 0.0;
	mm_tbuf_t *b = s->buf[tid];
	assert(s->n_seg[i] <= MM_MAX_SEG);
	if (mm_dbg_flag & MM_DBG_PRINT_QNAME) {
		fprintf(stderr, "QR\t%s\t%d\t%d\n", s->seq[off].name, tid, s->seq[off].l_seq);
		t = realtime();
	}
	for (j = 0; j < s->n_seg[i]; ++j) {
		if (s->n_seg[i] == 2 && ((j == 0 && (pe_ori>>1&1)) || (j == 1 && (pe_ori&1))))
			mm_revcomp_bseq(&s->seq[off + j]);
		qlens[j] = s->seq[off + j].l_seq;
		qseqs[j] = s->seq[off + j].seq;
	}
	if (s->p->opt->flag & MM_F_INDEPEND_SEG) {
		for (j = 0; j < s->n_seg[i]; ++j) {
			mm_map_frag(s->p->mi, 1, &qlens[j], &qseqs[j], &s->n_reg[off+j], &s->reg[off+j], b, s->p->opt, s->seq[off+j].name);
			s->rep_len[off + j] = b->rep_len;
			s->frag_gap[off + j] = b->frag_gap;
		}
	} else {
		mm_map_frag(s->p->mi, s->n_seg[i], qlens, qseqs, &s->n_reg[off], &s->reg[off], b, s->p->opt, s->seq[off].name);
		for (j = 0; j < s->n_seg[i]; ++j) {
			s->rep_len[off + j] = b->rep_len;
			s->frag_gap[off + j] = b->frag_gap;
		}
	}
	for (j = 0; j < s->n_seg[i]; ++j) // flip the query strand and coordinate to the original read strand
		if (s->n_seg[i] == 2 && ((j == 0 && (pe_ori>>1&1)) || (j == 1 && (pe_ori&1)))) {
			int k, t;
			mm_revcomp_bseq(&s->seq[off + j]);
			for (k = 0; k < s->n_reg[off + j]; ++k) {
				mm_reg1_t *r = &s->reg[off + j][k];
				t = r->qs;
				r->qs = qlens[j] - r->qe;
				r->qe = qlens[j] - t;
				r->rev = !r->rev;
				if (r->p) {
					if (r->p->trans_strand == 1) r->p->trans_strand = 2;
					else if (r->p->trans_strand == 2) r->p->trans_strand = 1;
				}
			}
		}
	if (mm_dbg_flag & MM_DBG_PRINT_QNAME)
		fprintf(stderr, "QT\t%s\t%d\t%.6f\n", s->seq[off].name, tid, realtime() - t);
}

static void merge_hits(step_t *s)
{
	int f, i, k0, k, max_seg = 0, *n_reg_part, *rep_len_part, *frag_gap_part, *qlens;
	void *km;
	FILE **fp = s->p->fp_parts;
	const mm_mapopt_t *opt = s->p->opt;

	km = km_init();
	for (f = 0; f < s->n_frag; ++f)
		max_seg = max_seg > s->n_seg[f]? max_seg : s->n_seg[f];
	qlens = CALLOC(int, max_seg + s->p->n_parts * 3);
	n_reg_part = qlens + max_seg;
	rep_len_part = n_reg_part + s->p->n_parts;
	frag_gap_part = rep_len_part + s->p->n_parts;
	for (f = 0, k = k0 = 0; f < s->n_frag; ++f) {
		k0 = k;
		for (i = 0; i < s->n_seg[f]; ++i, ++k) {
			int j, l, t, rep_len = 0;
			qlens[i] = s->seq[k].l_seq;
			for (j = 0, s->n_reg[k] = 0; j < s->p->n_parts; ++j) {
				mm_err_fread(&n_reg_part[j],	sizeof(int), 1, fp[j]);
				mm_err_fread(&rep_len_part[j],  sizeof(int), 1, fp[j]);
				mm_err_fread(&frag_gap_part[j], sizeof(int), 1, fp[j]);
				s->n_reg[k] += n_reg_part[j];
				if (rep_len < rep_len_part[j])
					rep_len = rep_len_part[j];
			}
			s->reg[k] = CALLOC(mm_reg1_t, s->n_reg[k]);
			for (j = 0, l = 0; j < s->p->n_parts; ++j) {
				for (t = 0; t < n_reg_part[j]; ++t, ++l) {
					mm_reg1_t *r = &s->reg[k][l];
					uint32_t capacity;
					mm_err_fread(r, sizeof(mm_reg1_t), 1, fp[j]);
					r->rid += s->p->rid_shift[j];
					if (opt->flag & MM_F_CIGAR) {
						mm_err_fread(&capacity, 4, 1, fp[j]);
						r->p = (mm_extra_t*)calloc(capacity, 4);
						r->p->capacity = capacity;
						mm_err_fread(r->p, r->p->capacity, 4, fp[j]);
					}
				}
			}
			if (!(opt->flag&MM_F_SR) && s->seq[k].l_seq >= opt->rank_min_len)
				mm_update_dp_max(s->seq[k].l_seq, s->n_reg[k], s->reg[k], opt->rank_frac, opt->a, opt->b);
			for (j = 0; j < s->n_reg[k]; ++j) {
				mm_reg1_t *r = &s->reg[k][j];
				if (r->p) r->p->dp_max2 = 0; // reset ->dp_max2 as mm_set_parent() doesn't clear it; necessary with mm_update_dp_max()
				r->subsc = 0; // this may not be necessary
				r->n_sub = 0; // n_sub will be an underestimate as we don't see all the chains now, but it can't be accurate anyway
			}
			mm_hit_sort(km, &s->n_reg[k], s->reg[k], opt->alt_drop);
			mm_set_parent(km, opt->mask_level, opt->mask_len, s->n_reg[k], s->reg[k], opt->a * 2 + opt->b, opt->flag&MM_F_HARD_MLEVEL, opt->alt_drop);
			if (!(opt->flag & MM_F_ALL_CHAINS)) {
				mm_select_sub(km, opt->pri_ratio, s->p->mi->k*2, opt->best_n, 0, opt->max_gap * 0.8, &s->n_reg[k], s->reg[k]);
				mm_set_sam_pri(s->n_reg[k], s->reg[k]);
			}
			//s->reg[k] = remove_second_suboptimal_alignment(s->reg[k], &(s->n_reg[k]), k);

			mm_set_mapq2(km, s->n_reg[k], s->reg[k], opt->min_chain_score, opt->a, rep_len, !!(opt->flag & (MM_F_SR|MM_F_SR_RNA)), !!(opt->flag & MM_F_SPLICE), NULL);
		}
		if (s->n_seg[f] == 2 && opt->pe_ori >= 0 && (opt->flag&MM_F_CIGAR))
			mm_pair(km, frag_gap_part[0], opt->pe_bonus, opt->a * 2 + opt->b, opt->a, qlens, &s->n_reg[k0], &s->reg[k0], !!(opt->flag & MM_F_MATE_RESOLVE));
	}
	free(qlens);
	km_destroy(km);
}

static void *worker_pipeline(void *shared, int step, void *in)
{
	int i, j, k;
	pipeline_t *p = (pipeline_t*)shared;
	if (step == 0) { // step 0: read sequences
		int with_qual = (!!(p->opt->flag & MM_F_OUT_SAM) && !(p->opt->flag & MM_F_NO_QUAL));
		int with_comment = !!(p->opt->flag & MM_F_COPY_COMMENT);
		int frag_mode = (p->n_fp > 1 || !!(p->opt->flag & MM_F_FRAG_MODE));
		step_t *s;
		s = (step_t*)calloc(1, sizeof(step_t));
		if (p->n_fp > 1) s->seq = mm_bseq_read_frag2(p->n_fp, p->fp, p->mini_batch_size, with_qual, with_comment, &s->n_seq);
		else s->seq = mm_bseq_read3(p->fp[0], p->mini_batch_size, with_qual, with_comment, frag_mode, &s->n_seq);
		if (s->seq) {
			s->p = p;
			for (i = 0; i < s->n_seq; ++i)
				s->seq[i].rid = p->n_processed++;
			s->buf = (mm_tbuf_t**)calloc(p->n_threads, sizeof(mm_tbuf_t*));
			for (i = 0; i < p->n_threads; ++i)
				s->buf[i] = mm_tbuf_init();
			s->n_reg = (int*)calloc(5 * s->n_seq, sizeof(int));
			s->seg_off = s->n_reg + s->n_seq; // seg_off, n_seg, rep_len and frag_gap are allocated together with n_reg
			s->n_seg = s->seg_off + s->n_seq;
			s->rep_len = s->n_seg + s->n_seq;
			s->frag_gap = s->rep_len + s->n_seq;
			s->reg = (mm_reg1_t**)calloc(s->n_seq, sizeof(mm_reg1_t*));
			for (i = 1, j = 0; i <= s->n_seq; ++i)
				if (i == s->n_seq || !frag_mode || !mm_qname_same(s->seq[i-1].name, s->seq[i].name)) {
					s->n_seg[s->n_frag] = i - j;
					s->seg_off[s->n_frag++] = j;
					j = i;
				}
			return s;
		} else free(s);
	} else if (step == 1) { // step 1: map
		if (p->n_parts > 0) merge_hits((step_t*)in);
		else kt_for(p->n_threads, worker_for, in, ((step_t*)in)->n_frag);
		return in;
	} else if (step == 2) { // step 2: output
		void *km = 0;
		step_t *s = (step_t*)in;
		const mm_idx_t *mi = p->mi;
		for (i = 0; i < p->n_threads; ++i) mm_tbuf_destroy(s->buf[i]);
		free(s->buf);
		if ((p->opt->flag & MM_F_OUT_CS) && !(mm_dbg_flag & MM_DBG_NO_KALLOC)) km = km_init();
		for (k = 0; k < s->n_frag; ++k) {
			int seg_st = s->seg_off[k], seg_en = s->seg_off[k] + s->n_seg[k];
			for (i = seg_st; i < seg_en; ++i) {
				mm_bseq1_t *t = &s->seq[i];
				if (p->opt->split_prefix && p->n_parts == 0) { // then write to temporary files
					mm_err_fwrite(&s->n_reg[i],    sizeof(int), 1, p->fp_split);
					mm_err_fwrite(&s->rep_len[i],  sizeof(int), 1, p->fp_split);
					mm_err_fwrite(&s->frag_gap[i], sizeof(int), 1, p->fp_split);
					for (j = 0; j < s->n_reg[i]; ++j) {
						mm_reg1_t *r = &s->reg[i][j];
						mm_err_fwrite(r, sizeof(mm_reg1_t), 1, p->fp_split);
						if (p->opt->flag & MM_F_CIGAR) {
							mm_err_fwrite(&r->p->capacity, 4, 1, p->fp_split);
							mm_err_fwrite(r->p, r->p->capacity, 4, p->fp_split);
						}
					}
				} else if (p->opt->flag & MM_F_OUT_JUNC) { // extra logic for --write-junc
					for (j = 0; j < s->n_reg[i]; ++j) {
						const mm_reg1_t *r = &s->reg[i][j];
						if (r->id != r->parent || r->mapq < 10) continue;
						mm_write_junc(&p->str, mi, t, r);
						if (p->str.l > 0) mm_err_puts(p->str.s);
					}
				} else if (s->n_reg[i] > 0) { // the query has at least one hit
					for (j = 0; j < s->n_reg[i]; ++j) {
						const mm_reg1_t *r = &s->reg[i][j];
						assert(!r->sam_pri || r->id == r->parent);
						if ((p->opt->flag & MM_F_NO_PRINT_2ND) && r->id != r->parent)
							continue;
						if (p->opt->flag & MM_F_OUT_SAM)
							mm_write_sam3(&p->str, mi, t, i - seg_st, j, s->n_seg[k], &s->n_reg[seg_st], (const mm_reg1_t*const*)&s->reg[seg_st], km, p->opt->flag, s->rep_len[i]);
						else
							mm_write_paf4(&p->str, mi, t, r, km, p->opt->flag, s->rep_len[i], s->n_seg[k], i - seg_st);
						mm_err_puts(p->str.s);
					}
				} else if ((p->opt->flag & MM_F_PAF_NO_HIT) || ((p->opt->flag & MM_F_OUT_SAM) && !(p->opt->flag & MM_F_SAM_HIT_ONLY))) { // output an empty hit, if requested
					if (p->opt->flag & MM_F_OUT_SAM)
						mm_write_sam3(&p->str, mi, t, i - seg_st, -1, s->n_seg[k], &s->n_reg[seg_st], (const mm_reg1_t*const*)&s->reg[seg_st], km, p->opt->flag, s->rep_len[i]);
					else
						mm_write_paf4(&p->str, mi, t, 0, 0, p->opt->flag, s->rep_len[i], s->n_seg[k], i - seg_st);
					mm_err_puts(p->str.s);
				}
			}
			for (i = seg_st; i < seg_en; ++i) {
				for (j = 0; j < s->n_reg[i]; ++j) free(s->reg[i][j].p);
				free(s->reg[i]);
				free(s->seq[i].seq); free(s->seq[i].name);
				if (s->seq[i].qual) free(s->seq[i].qual);
				if (s->seq[i].comment) free(s->seq[i].comment);
			}
		}
		free(s->reg); free(s->n_reg); free(s->seq); // seg_off, n_seg, rep_len and frag_gap were allocated with reg; no memory leak here
		km_destroy(km);
		if (mm_verbose >= 3)
			fprintf(stderr, "[M::%s::%.3f*%.2f] mapped %d sequences\n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), s->n_seq);
		free(s);
	}
	return 0;
}

static mm_bseq_file_t **open_bseqs(int n, const char **fn)
{
	mm_bseq_file_t **fp;
	int i, j;
	fp = (mm_bseq_file_t**)calloc(n, sizeof(mm_bseq_file_t*));
	for (i = 0; i < n; ++i) {
		if ((fp[i] = mm_bseq_open(fn[i])) == 0) {
			if (mm_verbose >= 1)
				fprintf(stderr, "ERROR: failed to open file '%s': %s\n", fn[i], strerror(errno));
			for (j = 0; j < i; ++j)
				mm_bseq_close(fp[j]);
			free(fp);
			return 0;
		}
	}
	return fp;
}

int mm_map_file_frag(const mm_idx_t *idx, int n_segs, const char **fn, const mm_mapopt_t *opt, int n_threads)
{
	int i, pl_threads;
	pipeline_t pl;
	if (n_segs < 1) return -1;
	memset(&pl, 0, sizeof(pipeline_t));
	pl.n_fp = n_segs;
	pl.fp = open_bseqs(pl.n_fp, fn);
	if (pl.fp == 0) return -1;
	pl.opt = opt, pl.mi = idx;
	pl.n_threads = n_threads > 1? n_threads : 1;
	pl.mini_batch_size = opt->mini_batch_size;
	if (opt->split_prefix)
		pl.fp_split = mm_split_init(opt->split_prefix, idx);
	pl_threads = n_threads == 1? 1 : (opt->flag&MM_F_2_IO_THREADS)? 3 : 2;
	kt_pipeline(pl_threads, worker_pipeline, &pl, 3);

	free(pl.str.s);
	if (pl.fp_split) fclose(pl.fp_split);
	for (i = 0; i < pl.n_fp; ++i)
		mm_bseq_close(pl.fp[i]);
	free(pl.fp);
	return 0;
}

int mm_map_file(const mm_idx_t *idx, const char *fn, const mm_mapopt_t *opt, int n_threads)
{
	return mm_map_file_frag(idx, 1, &fn, opt, n_threads);
}

int mm_split_merge(int n_segs, const char **fn, const mm_mapopt_t *opt, int n_split_idx)
{
	int i;
	pipeline_t pl;
	mm_idx_t *mi;
	if (n_segs < 1 || n_split_idx < 1) return -1;
	memset(&pl, 0, sizeof(pipeline_t));
	pl.n_fp = n_segs;
	pl.fp = open_bseqs(pl.n_fp, fn);
	if (pl.fp == 0) return -1;
	pl.opt = opt;
	pl.mini_batch_size = opt->mini_batch_size;

	pl.n_parts = n_split_idx;
	pl.fp_parts  = CALLOC(FILE*, pl.n_parts);
	pl.rid_shift = CALLOC(uint32_t, pl.n_parts);
	pl.mi = mi = mm_split_merge_prep(opt->split_prefix, n_split_idx, pl.fp_parts, pl.rid_shift);
	if (pl.mi == 0) {
		free(pl.fp_parts);
		free(pl.rid_shift);
		return -1;
	}
	for (i = n_split_idx - 1; i > 0; --i)
		pl.rid_shift[i] = pl.rid_shift[i - 1];
	for (pl.rid_shift[0] = 0, i = 1; i < n_split_idx; ++i)
		pl.rid_shift[i] += pl.rid_shift[i - 1];
	if (opt->flag & MM_F_OUT_SAM)
		for (i = 0; i < (int32_t)pl.mi->n_seq; ++i)
			printf("@SQ\tSN:%s\tLN:%d\n", pl.mi->seq[i].name, pl.mi->seq[i].len);

	kt_pipeline(2, worker_pipeline, &pl, 3);

	free(pl.str.s);
	mm_idx_destroy(mi);
	free(pl.rid_shift);
	for (i = 0; i < n_split_idx; ++i)
		fclose(pl.fp_parts[i]);
	free(pl.fp_parts);
	for (i = 0; i < pl.n_fp; ++i)
		mm_bseq_close(pl.fp[i]);
	free(pl.fp);
	mm_split_rm_tmp(opt->split_prefix, n_split_idx);
	return 0;
}
