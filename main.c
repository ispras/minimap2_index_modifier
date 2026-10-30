#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "bseq.h"
#include "minimap.h"
#include "mmpriv.h"
#include "mm_lift.h"
#include "mm_varalt.h"
#include "mm_hapdp.h"
#include "rmin.h"
#include "ketopt.h"

#ifdef __linux__
#include <sys/resource.h>
#include <sys/time.h>
void liftrlimit()
{
	struct rlimit r;
	getrlimit(RLIMIT_AS, &r);
	r.rlim_cur = r.rlim_max;
	setrlimit(RLIMIT_AS, &r);
}
#else
void liftrlimit() {}
#endif

static ko_longopt_t long_options[] = {
	{ "bucket-bits",    ko_required_argument, 300 },
	{ "mb-size",        ko_required_argument, 'K' },
	{ "seed",           ko_required_argument, 302 },
	{ "no-kalloc",      ko_no_argument,       303 },
	{ "print-qname",    ko_no_argument,       304 },
	{ "no-self",        ko_no_argument,       'D' },
	{ "print-seeds",    ko_no_argument,       306 },
	{ "max-chain-skip", ko_required_argument, 307 },
	{ "min-dp-len",     ko_required_argument, 308 },
	{ "print-aln-seq",  ko_no_argument,       309 },
	{ "splice",         ko_no_argument,       310 },
	{ "cost-non-gt-ag", ko_required_argument, 'C' },
	{ "no-long-join",   ko_no_argument,       312 },
	{ "sr",             ko_optional_argument, 313 },
	{ "frag",           ko_required_argument, 314 },
	{ "secondary",      ko_required_argument, 315 },
	{ "cs",             ko_optional_argument, 316 },
	{ "end-bonus",      ko_required_argument, 317 },
	{ "no-pairing",     ko_no_argument,       318 }, // deprecated but reserved for backward compatibility
	{ "splice-flank",   ko_required_argument, 319 },
	{ "idx-no-seq",     ko_no_argument,       320 },
	{ "end-seed-pen",   ko_required_argument, 321 },
	{ "for-only",       ko_no_argument,       322 },
	{ "rev-only",       ko_no_argument,       323 },
	{ "heap-sort",      ko_required_argument, 324 },
	{ "all-chain",      ko_no_argument,       'P' },
	{ "dual",           ko_required_argument, 326 },
	{ "max-clip-ratio", ko_required_argument, 327 },
	{ "min-occ-floor",  ko_required_argument, 328 },
	{ "MD",             ko_no_argument,       329 },
	{ "lj-min-ratio",   ko_required_argument, 330 },
	{ "score-N",        ko_required_argument, 331 },
	{ "eqx",            ko_no_argument,       332 },
	{ "paf-no-hit",     ko_no_argument,       333 },
	{ "split-prefix",   ko_required_argument, 334 },
	{ "no-end-flt",     ko_no_argument,       335 },
	{ "hard-mask-level",ko_no_argument,       336 },
	{ "cap-sw-mem",     ko_required_argument, 337 },
	{ "max-qlen",       ko_required_argument, 338 },
	{ "max-chain-iter", ko_required_argument, 339 },
	{ "junc-bed",       ko_required_argument, 340 },
	{ "junc-bonus",     ko_required_argument, 341 },
	{ "sam-hit-only",   ko_no_argument,       342 },
	{ "chain-gap-scale",ko_required_argument, 343 },
	{ "alt",            ko_required_argument, 344 },
	{ "alt-drop",       ko_required_argument, 345 },
	{ "mask-len",       ko_required_argument, 346 },
	{ "rmq",            ko_optional_argument, 347 },
	{ "qstrand",        ko_no_argument,       348 },
	{ "cap-kalloc",     ko_required_argument, 349 },
	{ "q-occ-frac",     ko_required_argument, 350 },
	{ "chain-skip-scale",ko_required_argument,351 },
	{ "print-chains",   ko_no_argument,       352 },
	{ "no-hash-name",   ko_no_argument,       353 },
	{ "secondary-seq",  ko_no_argument,       354 },
	{ "ds",             ko_no_argument,       355 },
	{ "rmq-inner",      ko_required_argument, 356 },
	{ "spsc",           ko_required_argument, 357 },
	{ "junc-pen",       ko_required_argument, 358 },
	{ "pairing",        ko_required_argument, 359 },
	{ "jump-min-match", ko_required_argument, 360 },
	{ "write-junc",     ko_no_argument,       361 },
	{ "pass1",          ko_required_argument, 362 },
	{ "spsc-scale",     ko_required_argument, 363 },
	{ "spsc0",          ko_required_argument, 364 },
	{ "dbg-seed-occ",   ko_no_argument,       501 },
    { "vcf-file-with-variants",ko_required_argument, 700 },
	{ "mapq-quotient", ko_no_argument,       710 },
	{ "mapq-prefer-primary", ko_no_argument, 711 },
	{ "lift",           ko_no_argument,       712 },
	{ "lift-map",       ko_required_argument, 713 },
	{ "lift-junction",  ko_required_argument, 714 },
	{ "mate-resolve",   ko_no_argument,       715 },
	{ "mate-mapq",      ko_required_argument, 716 },
	{ "lift-noshift",   ko_no_argument,       717 },
	{ "no-mapq-rule",   ko_no_argument,       718 },
	{ "lift-realign",   ko_no_argument,       719 },
	{ "lift-tie-base",  ko_no_argument,       720 },
	{ "lift-tie",       ko_required_argument, 721 },
	{ "mate-frag",      ko_required_argument, 722 },
	{ "pair-mapq",      ko_required_argument, 770 },
	{ "mate-sd",        ko_required_argument, 771 },
	{ "mate-gap",       ko_required_argument, 772 },
	{ "lift-realign-margin", ko_required_argument, 723 },
	{ "var-alt",        ko_required_argument, 724 },
	{ "hap-dp",         ko_required_argument, 760 },
	{ "hap-dp-max",     ko_required_argument, 761 },
	{ "hap-dp-mask",    ko_required_argument, 762 },
    { "parse-haplotype",ko_no_argument,       701 },
	{ "read-minimizers", ko_required_argument, 730 },
	{ "read-kmers",      ko_required_argument, 731 },
	{ "read-min-count",  ko_required_argument, 732 },
	{ "read-min-mode",   ko_required_argument, 733 },
	{ "read-min-frac",   ko_required_argument, 734 },
	{ "read-min-rate",   ko_required_argument, 735 },
	{ "read-min-seed",   ko_required_argument, 736 },
	{ "read-min-list",   ko_required_argument, 737 },
	{ "read-min-dump",   ko_required_argument, 738 },
	{ "read-min-ref-kmers", ko_required_argument, 739 },
	{ "help",           ko_no_argument,       'h' },
	{ "max-intron-len", ko_required_argument, 'G' },
	{ "version",        ko_no_argument,       'V' },
	{ "min-count",      ko_required_argument, 'n' },
	{ "min-chain-score",ko_required_argument, 'm' },
	{ "mask-level",     ko_required_argument, 'M' },
	{ "min-dp-score",   ko_required_argument, 's' },
	{ "sam",            ko_no_argument,       'a' },
	{ 0, 0, 0 }
};

static inline int64_t mm_parse_num2(const char *str, char **q)
{
	double x;
	char *p;
	x = strtod(str, &p);
	if (*p == 'G' || *p == 'g') x *= 1e9, ++p;
	else if (*p == 'M' || *p == 'm') x *= 1e6, ++p;
	else if (*p == 'K' || *p == 'k') x *= 1e3, ++p;
	if (q) *q = p;
	return (int64_t)(x + .499);
}

static inline int64_t mm_parse_num(const char *str)
{
	return mm_parse_num2(str, 0);
}

static inline void yes_or_no(mm_mapopt_t *opt, int64_t flag, int long_idx, const char *arg, int yes_to_set)
{
	if (yes_to_set) {
		if (strcmp(arg, "yes") == 0 || strcmp(arg, "y") == 0) opt->flag |= flag;
		else if (strcmp(arg, "no") == 0 || strcmp(arg, "n") == 0) opt->flag &= ~flag;
		else fprintf(stderr, "[WARNING]\033[1;31m option '--%s' only accepts 'yes' or 'no'.\033[0m\n", long_options[long_idx].name);
	} else {
		if (strcmp(arg, "yes") == 0 || strcmp(arg, "y") == 0) opt->flag &= ~flag;
		else if (strcmp(arg, "no") == 0 || strcmp(arg, "n") == 0) opt->flag |= flag;
		else fprintf(stderr, "[WARNING]\033[1;31m option '--%s' only accepts 'yes' or 'no'.\033[0m\n", long_options[long_idx].name);
	}
}

int main(int argc, char *argv[])
{
	const char *opt_str = "2aSDw:k:K:t:r:f:Vv:g:G:I:d:XT:s:x:Hcp:M:n:z:A:B:b:O:E:m:N:Qu:R:hF:LC:yYPo:e:U:J:j:";
	ketopt_t o = KETOPT_INIT;
	mm_mapopt_t opt;
	mm_idxopt_t ipt;
	int i, c, n_threads = 3, n_parts, old_best_n = -1;
	char *vcf_with_variants = 0;
	char *fn_lift_map = 0;
	char *fn_var_alt = 0, *fn_hap_dp = 0, *fn_hap_dp_mask = 0;
	int hap_dp_max = 16;
	int lift_jct = 101;
	float spsc_scale = 0.7f;
	char *fnw = 0, *rg = 0, *fn_bed_junc = 0, *fn_bed_jump = 0, *fn_bed_pass1 = 0, *fn_spsc = 0, *s, *alt_list = 0;
	FILE *fp_help = stderr;
	mm_idx_reader_t *idx_rdr;
	mm_idx_t *mi;

	mm_verbose = 3;
	liftrlimit();
	mm_realtime0 = realtime();
	mm_set_opt(0, &ipt, &opt);

	while ((c = ketopt(&o, argc, argv, 1, opt_str, long_options)) >= 0) { // test command line options and apply option -x/preset first
		if (c == 'x') {
			if (mm_set_opt(o.arg, &ipt, &opt) < 0) {
				fprintf(stderr, "[ERROR] unknown preset '%s'\n", o.arg);
				return 1;
			}
		} else if (c == ':') {
			fprintf(stderr, "[ERROR] missing option argument\n");
			return 1;
		} else if (c == '?') {
			fprintf(stderr, "[ERROR] unknown option in \"%s\"\n", argv[o.i - 1]);
			return 1;
		}
	}
	o = KETOPT_INIT;

	while ((c = ketopt(&o, argc, argv, 1, opt_str, long_options)) >= 0) {
		if (c == 'w') ipt.w = atoi(o.arg);
		else if (c == 'k') ipt.k = atoi(o.arg);
		else if (c == 'H') ipt.flag |= MM_I_HPC;
		else if (c == 'd') fnw = o.arg; // the above are indexing related options, except -I
		else if (c == 't') n_threads = atoi(o.arg);
		else if (c == 'v') mm_verbose = atoi(o.arg);
		else if (c == 'g') opt.max_gap = (int)mm_parse_num(o.arg);
		else if (c == 'G') mm_mapopt_max_intron_len(&opt, (int)mm_parse_num(o.arg));
		else if (c == 'F') opt.max_frag_len = (int)mm_parse_num(o.arg);
		else if (c == 'N') old_best_n = opt.best_n, opt.best_n = atoi(o.arg);
		else if (c == 'p') opt.pri_ratio = atof(o.arg);
		else if (c == 'M') opt.mask_level = atof(o.arg);
		else if (c == 'c') opt.flag |= MM_F_OUT_CG | MM_F_CIGAR;
		else if (c == 'D') opt.flag |= MM_F_NO_DIAG;
		else if (c == 'P') opt.flag |= MM_F_ALL_CHAINS;
		else if (c == 'X') opt.flag |= MM_F_ALL_CHAINS | MM_F_NO_DIAG | MM_F_NO_DUAL | MM_F_NO_LJOIN; // -D -P --no-long-join --dual=no
		else if (c == 'a') opt.flag |= MM_F_OUT_SAM | MM_F_CIGAR;
		else if (c == 'Q') opt.flag |= MM_F_NO_QUAL;
		else if (c == 'Y') opt.flag |= MM_F_SOFTCLIP;
		else if (c == 'L') opt.flag |= MM_F_LONG_CIGAR;
		else if (c == 'y') opt.flag |= MM_F_COPY_COMMENT;
		else if (c == 'T') opt.sdust_thres = atoi(o.arg);
		else if (c == 'n') opt.min_cnt = atoi(o.arg);
		else if (c == 'm') opt.min_chain_score = atoi(o.arg);
		else if (c == 'A') opt.a = atoi(o.arg);
		else if (c == 'B') opt.b = atoi(o.arg);
		else if (c == 'b') opt.transition = atoi(o.arg);
		else if (c == 's') opt.min_dp_max = atoi(o.arg);
		else if (c == 'C') opt.noncan = atoi(o.arg);
		else if (c == 'I') ipt.batch_size = mm_parse_num(o.arg);
		else if (c == 'K') opt.mini_batch_size = mm_parse_num(o.arg);
		else if (c == 'e') opt.occ_dist = mm_parse_num(o.arg);
		else if (c == 'R') rg = o.arg;
		else if (c == 'h') fp_help = stdout;
		else if (c == '2') opt.flag |= MM_F_2_IO_THREADS;
		else if (c == 'j') fn_bed_jump = o.arg;
		else if (c == 'J') {
			int t;
			t = atoi(o.arg);
			if (t == 0) opt.flag |= MM_F_SPLICE_OLD;
			else if (t == 1) opt.flag &= ~MM_F_SPLICE_OLD;
		} else if (c == 'o') {
			if (strcmp(o.arg, "-") != 0) {
				if (freopen(o.arg, "wb", stdout) == NULL) {
					fprintf(stderr, "[ERROR]\033[1;31m failed to write the output to file '%s'\033[0m: %s\n", o.arg, strerror(errno));
					exit(1);
				}
			}
		}
		else if (c == 300) ipt.bucket_bits = atoi(o.arg); // --bucket-bits
		else if (c == 302) opt.seed = atoi(o.arg); // --seed
		else if (c == 303) mm_dbg_flag |= MM_DBG_NO_KALLOC; // --no-kalloc
		else if (c == 304) mm_dbg_flag |= MM_DBG_PRINT_QNAME; // --print-qname
		else if (c == 306) mm_dbg_flag |= MM_DBG_PRINT_QNAME | MM_DBG_PRINT_SEED, n_threads = 1; // --print-seed
		else if (c == 307) opt.max_chain_skip = atoi(o.arg); // --max-chain-skip
		else if (c == 339) opt.max_chain_iter = atoi(o.arg); // --max-chain-iter
		else if (c == 308) opt.min_ksw_len = atoi(o.arg); // --min-dp-len
		else if (c == 309) mm_dbg_flag |= MM_DBG_PRINT_QNAME | MM_DBG_PRINT_ALN_SEQ, n_threads = 1; // --print-aln-seq
		else if (c == 310) opt.flag |= MM_F_SPLICE; // --splice
		else if (c == 312) opt.flag |= MM_F_NO_LJOIN; // --no-long-join
		else if (c == 317) opt.end_bonus = atoi(o.arg); // --end-bonus
		else if (c == 318) opt.flag |= MM_F_INDEPEND_SEG; // --no-pairing (deprecated)
		else if (c == 320) ipt.flag |= MM_I_NO_SEQ; // --idx-no-seq
		else if (c == 321) opt.anchor_ext_shift = atoi(o.arg); // --end-seed-pen
		else if (c == 322) opt.flag |= MM_F_FOR_ONLY; // --for-only
		else if (c == 323) opt.flag |= MM_F_REV_ONLY; // --rev-only
		else if (c == 327) opt.max_clip_ratio = atof(o.arg); // --max-clip-ratio
		else if (c == 328) opt.min_mid_occ = atoi(o.arg); // --min-occ-floor
		else if (c == 329) opt.flag |= MM_F_OUT_MD; // --MD
		else if (c == 331) opt.sc_ambi = atoi(o.arg); // --score-N
		else if (c == 332) opt.flag |= MM_F_EQX; // --eqx
		else if (c == 333) opt.flag |= MM_F_PAF_NO_HIT; // --paf-no-hit
		else if (c == 334) opt.split_prefix = o.arg; // --split-prefix
		else if (c == 335) opt.flag |= MM_F_NO_END_FLT; // --no-end-flt
		else if (c == 336) opt.flag |= MM_F_HARD_MLEVEL; // --hard-mask-level
		else if (c == 337) opt.max_sw_mat = mm_parse_num(o.arg); // --cap-sw-mat
		else if (c == 338) opt.max_qlen = mm_parse_num(o.arg); // --max-qlen
		else if (c == 340) fn_bed_junc = o.arg; // --junc-bed
		else if (c == 341) opt.junc_bonus = atoi(o.arg); // --junc-bonus
		else if (c == 342) opt.flag |= MM_F_SAM_HIT_ONLY; // --sam-hit-only
		else if (c == 343) opt.chain_gap_scale = atof(o.arg); // --chain-gap-scale
		else if (c == 351) opt.chain_skip_scale = atof(o.arg); // --chain-skip-scale
		else if (c == 344) alt_list = o.arg; // --alt
		else if (c == 345) opt.alt_drop = atof(o.arg); // --alt-drop
		else if (c == 346) opt.mask_len = mm_parse_num(o.arg); // --mask-len
		else if (c == 348) opt.flag |= MM_F_QSTRAND | MM_F_NO_INV; // --qstrand
		else if (c == 349) opt.cap_kalloc = mm_parse_num(o.arg); // --cap-kalloc
		else if (c == 350) opt.q_occ_frac = atof(o.arg); // --q-occ-frac
		else if (c == 352) mm_dbg_flag |= MM_DBG_PRINT_CHAIN; // --print-chains
		else if (c == 353) opt.flag |= MM_F_NO_HASH_NAME; // --no-hash-name
		else if (c == 354) opt.flag |= MM_F_SECONDARY_SEQ; // --secondary-seq
		else if (c == 355) opt.flag |= MM_F_OUT_DS; // --ds
		else if (c == 356) opt.rmq_inner_dist = mm_parse_num(o.arg); // --rmq-inner
		else if (c == 357) fn_spsc = o.arg; // --spsc
		else if (c == 360) opt.jump_min_match = mm_parse_num(o.arg); // --jump-min-match
		else if (c == 361) opt.flag |= MM_F_OUT_JUNC | MM_F_CIGAR; // --write-junc
		else if (c == 362) fn_bed_pass1 = o.arg; // --jump-pass1
		else if (c == 501) mm_dbg_flag |= MM_DBG_SEED_FREQ; // --dbg-seed-occ
		else if (c == 363) spsc_scale = atof(o.arg); // --spsc-scale
		else if (c == 358 || c == 364) opt.junc_pen = atoi(o.arg); // --junc-pen or --spsc0
        else if (c == 700) vcf_with_variants = o.arg; // --vcf-file-with-variants
		else if (c == 710) opt.flag |= MM_F_MAPQ_QUOTIENT; // --mapq-quotient
		else if (c == 711) opt.flag |= MM_F_MAPQ_PREFER_T; // --mapq-prefer-primary
		else if (c == 712) opt.flag |= MM_F_LIFT; // --lift
		else if (c == 713) fn_lift_map = o.arg; // --lift-map
		else if (c == 714) lift_jct = atoi(o.arg); // --lift-junction
		else if (c == 715) opt.flag |= MM_F_MATE_RESOLVE; // --mate-resolve
		else if (c == 716) opt.mate_mapq = atoi(o.arg); // --mate-mapq
		else if (c == 717) opt.flag |= MM_F_LIFT_NOSHIFT; // --lift-noshift
		else if (c == 718) opt.flag |= MM_F_NO_MAPQ_RULE; // --no-mapq-rule
		else if (c == 719) opt.flag |= MM_F_LIFT_REALIGN; // --lift-realign
		else if (c == 720) opt.flag |= MM_F_LIFT_TIE_BASE; // --lift-tie-base
		else if (c == 721) opt.lift_tie = atoi(o.arg); // --lift-tie
		else if (c == 722) opt.mate_frag = atoi(o.arg); // --mate-frag
		else if (c == 770) { // --pair-mapq
			if (strcmp(o.arg, "off") == 0) opt.pair_mapq = 0;
			else if (strcmp(o.arg, "vg") == 0) opt.pair_mapq = 1;
			else if (strcmp(o.arg, "vgup") == 0) opt.pair_mapq = 2;
			else { fprintf(stderr, "[ERROR] --pair-mapq takes off, vg or vgup\n"); return 1; }
		}
		else if (c == 771) opt.mate_sd = atoi(o.arg); // --mate-sd
		else if (c == 772) opt.mate_gap = atoi(o.arg); // --mate-gap
		else if (c == 723) opt.lift_realign_margin = atoi(o.arg); // --lift-realign-margin
		else if (c == 724) fn_var_alt = o.arg; // --var-alt
		else if (c == 760) fn_hap_dp = o.arg; // --hap-dp
		else if (c == 761) hap_dp_max = atoi(o.arg); // --hap-dp-max
		else if (c == 762) fn_hap_dp_mask = o.arg; // --hap-dp-mask
        else if (c == 701) ipt.flag |= MM_PARSE_HT;  // --parse-haplotype
		else if (c == 730) rmin_opt.fn_mm = o.arg; // --read-minimizers (rmin)
		else if (c == 731) rmin_opt.fn_km = o.arg; // --read-kmers
		else if (c == 732) rmin_opt.min_cnt = atoi(o.arg); // --read-min-count
		else if (c == 733) { // --read-min-mode
			if ((rmin_opt.mode = rmin_parse_mode(o.arg)) < 0) { fprintf(stderr, "[ERROR] unknown --read-min-mode '%s'\n", o.arg); return 1; }
		}
		else if (c == 734) rmin_opt.frac = atof(o.arg); // --read-min-frac
		else if (c == 735) rmin_opt.rate = atof(o.arg); // --read-min-rate
		else if (c == 736) rmin_opt.seed = strtoull(o.arg, 0, 10); // --read-min-seed
		else if (c == 737) rmin_opt.fn_list = o.arg; // --read-min-list
		else if (c == 738) rmin_opt.fn_dump = o.arg; // --read-min-dump
		else if (c == 739) rmin_opt.fn_ref = o.arg; // --read-min-ref-kmers
		else if (c == 330) {
			fprintf(stderr, "[WARNING] \033[1;31m --lj-min-ratio has been deprecated.\033[0m\n");
		} else if (c == 313) { // --sr
			if (o.arg == 0 || strcmp(o.arg, "dna") == 0) {
				opt.flag |= MM_F_SR;
			} else if (strcmp(o.arg, "rna") == 0) {
				opt.flag |= MM_F_SR_RNA;
			} else if (strcmp(o.arg, "no") == 0) {
				opt.flag &= ~(uint64_t)(MM_F_SR|MM_F_SR_RNA);
			} else if (mm_verbose >= 2) {
				opt.flag |= MM_F_SR;
				fprintf(stderr, "[WARNING]\033[1;31m --sr only takes 'dna' or 'rna'. Invalid values are assumed to be 'dna'.\033[0m\n");
			}
		} else if (c == 314) { // --frag
			yes_or_no(&opt, MM_F_FRAG_MODE, o.longidx, o.arg, 1);
		} else if (c == 315) { // --secondary
			yes_or_no(&opt, MM_F_NO_PRINT_2ND, o.longidx, o.arg, 0);
		} else if (c == 316) { // --cs
			opt.flag |= MM_F_OUT_CS | MM_F_CIGAR;
			if (o.arg == 0 || strcmp(o.arg, "short") == 0) {
				opt.flag &= ~MM_F_OUT_CS_LONG;
			} else if (strcmp(o.arg, "long") == 0) {
				opt.flag |= MM_F_OUT_CS_LONG;
			} else if (strcmp(o.arg, "none") == 0) {
				opt.flag &= ~MM_F_OUT_CS;
			} else if (mm_verbose >= 2) {
				fprintf(stderr, "[WARNING]\033[1;31m --cs only takes 'short' or 'long'. Invalid values are assumed to be 'short'.\033[0m\n");
			}
		} else if (c == 319) { // --splice-flank
			yes_or_no(&opt, MM_F_SPLICE_FLANK, o.longidx, o.arg, 1);
		} else if (c == 324) { // --heap-sort
			yes_or_no(&opt, MM_F_HEAP_SORT, o.longidx, o.arg, 1);
		} else if (c == 326) { // --dual
			yes_or_no(&opt, MM_F_NO_DUAL, o.longidx, o.arg, 0);
		} else if (c == 347) { // --rmq
			if (o.arg) yes_or_no(&opt, MM_F_RMQ, o.longidx, o.arg, 1);
			else opt.flag |= MM_F_RMQ;
		} else if (c == 359) { // --pairing
			if (strcmp(o.arg, "no") == 0) opt.flag |= MM_F_INDEPEND_SEG;
			else if (strcmp(o.arg, "weak") == 0) opt.flag |= MM_F_WEAK_PAIRING, opt.flag &= ~(uint64_t)MM_F_INDEPEND_SEG;
			else {
				if (strcmp(o.arg, "strong") != 0 && mm_verbose >= 2)
					fprintf(stderr, "[WARNING]\033[1;31m unrecognized argument for --pairing; assuming 'strong'.\033[0m\n");
				opt.flag &= ~(uint64_t)(MM_F_INDEPEND_SEG|MM_F_WEAK_PAIRING);
			}
		} else if (c == 'S') {
			opt.flag |= MM_F_OUT_CS | MM_F_CIGAR | MM_F_OUT_CS_LONG;
			if (mm_verbose >= 2)
				fprintf(stderr, "[WARNING]\033[1;31m option -S is deprecated and may be removed in future. Please use --cs=long instead.\033[0m\n");
		} else if (c == 'V') {
			puts(MM_VERSION);
			return 0;
		} else if (c == 'r') {
			opt.bw = (int)mm_parse_num2(o.arg, &s);
			if (*s == ',') opt.bw_long = (int)mm_parse_num2(s + 1, &s);
		} else if (c == 'U') {
			opt.min_mid_occ = strtol(o.arg, &s, 10);
			if (*s == ',') opt.max_mid_occ = strtol(s + 1, &s, 10);
		} else if (c == 'f') {
			double x;
			char *p;
			x = strtod(o.arg, &p);
			if (x < 1.0) opt.mid_occ_frac = x, opt.mid_occ = 0;
			else opt.mid_occ = (int)(x + .499);
			if (*p == ',') opt.max_occ = (int)(strtod(p+1, &p) + .499);
		} else if (c == 'u') {
			if (*o.arg == 'b') opt.flag |= MM_F_SPLICE_FOR|MM_F_SPLICE_REV; // both strands
			else if (*o.arg == 'f') opt.flag |= MM_F_SPLICE_FOR, opt.flag &= ~MM_F_SPLICE_REV; // match GT-AG
			else if (*o.arg == 'r') opt.flag |= MM_F_SPLICE_REV, opt.flag &= ~MM_F_SPLICE_FOR; // match CT-AC (reverse complement of GT-AG)
			else if (*o.arg == 'n') opt.flag &= ~(MM_F_SPLICE_FOR|MM_F_SPLICE_REV); // don't try to match the GT-AG signal
			else {
				fprintf(stderr, "[ERROR]\033[1;31m unrecognized cDNA direction\033[0m\n");
				return 1;
			}
		} else if (c == 'z') {
			opt.zdrop = opt.zdrop_inv = strtol(o.arg, &s, 10);
			if (*s == ',') opt.zdrop_inv = strtol(s + 1, &s, 10);
		} else if (c == 'O') {
			opt.q = opt.q2 = strtol(o.arg, &s, 10);
			if (*s == ',') opt.q2 = strtol(s + 1, &s, 10);
		} else if (c == 'E') {
			opt.e = opt.e2 = strtol(o.arg, &s, 10);
			if (*s == ',') opt.e2 = strtol(s + 1, &s, 10);
		}
	}
	if (!fnw && !(opt.flag&MM_F_CIGAR))
		ipt.flag |= MM_I_NO_SEQ;
	if (mm_check_opt(&ipt, &opt) < 0)
		return 1;
	if (opt.best_n == 0) {
		fprintf(stderr, "[WARNING]\033[1;31m changed '-N 0' to '-N %d --secondary=no'.\033[0m\n", old_best_n);
		opt.best_n = old_best_n, opt.flag |= MM_F_NO_PRINT_2ND;
	}

	if (argc == o.ind || fp_help == stdout) {
		fprintf(fp_help, "Usage: minimap2 [options] <target.fa>|<target.idx> [query.fa] [...]\n");
		fprintf(fp_help, "Options:\n");
		fprintf(fp_help, "  Indexing:\n");
		fprintf(fp_help, "    -H           use homopolymer-compressed k-mer (preferrable for PacBio)\n");
		fprintf(fp_help, "    -k INT       k-mer size (no larger than 28) [%d]\n", ipt.k);
		fprintf(fp_help, "    -w INT       minimizer window size [%d]\n", ipt.w);
		fprintf(fp_help, "    -I NUM       split index for every ~NUM input bases [8G]\n");
		fprintf(fp_help, "    -d FILE      dump index to FILE []\n");
		fprintf(fp_help, "    --vcf-file-with-variants FILE      pass VCF FILE to modify index []\n");
		fprintf(fp_help, "    --lift             emit base-reference coordinates for augmentation-contig\n");
		fprintf(fp_help, "                       placements (no external liftover)\n");
		fprintf(fp_help, "    --lift-map FILE    sidecar TSV for --lift: <contig>\\t<ins_len>\\t<del_len>[\\t<jct>] (or <contig>\\t<|ALT|-|REF|>)\n");
		fprintf(fp_help, "    --lift-junction INT fallback junction (left flank + 1) for sidecar lines without the 4th column [101]\n");
		fprintf(fp_help, "    --lift-noshift     (rollback only) no indel left-shift / NM refresh after --lift\n");
		fprintf(fp_help, "    --no-mapq-rule     (control only) no candidate-collapsing rule at all, as in the production image\n");
		fprintf(fp_help, "    --lift-realign     after --lift, re-derive each lifted read's CIGAR against the base as mm_align1 does:\n");
		fprintf(fp_help, "                       exact anchor from the projection + extension DP to both read ends (fitting where\n");
		fprintf(fp_help, "                       the projection reached the end); replaces the projection iff not worse by score\n");
		fprintf(fp_help, "    --lift-realign-margin INT  a re-derived middle replaces the projection only if better by MORE than\n");
		fprintf(fp_help, "                       this score (keeps one representation per locus) [a*2+b]; 0 = any improvement\n");
		fprintf(fp_help, "    --lift-tie-base    use a contig only where it is strictly better: if within --lift-tie of the best\n");
		fprintf(fp_help, "                       score there are two loci and one is a contig, drop the contig candidates and\n");
		fprintf(fp_help, "                       place the read as on the plain reference (tag tb:i:1)\n");
		fprintf(fp_help, "    --lift-tie INT     score window that counts as a tie for --lift-tie-base [a*2+b]\n");
		fprintf(fp_help, "    --var-alt FILE     known SNV alleles <seq>\\t<1-based pos>\\t<ALTs> (text or gzip): a mismatch whose read\n");
		fprintf(fp_help, "                       base is a known ALT at that base position is scored as a match (dp, identity for\n");
		fprintf(fp_help, "                       mapq; CIGAR and NM unchanged; tag va:i); candidates re-ranked before mapq\n");
		fprintf(fp_help, "    --hap-dp FILE      haplotype-consistent DP: <seq>\\t<1-based pos>\\t<ALT>\\t<mask_hi>\\t<mask_lo> (hex carrier\n");
		fprintf(fp_help, "                       haplotypes); each candidate is scored on the best panel haplotype string over\n");
		fprintf(fp_help, "                       the SNVs of its window, the base string included (dp, identity for mapq;\n");
		fprintf(fp_help, "                       CIGAR and NM against the base; tag hd:i); not with --var-alt\n");
		fprintf(fp_help, "    --hap-dp-max INT   haplotype strings per window, the base included, most frequent first [16]\n");
		fprintf(fp_help, "    --hap-dp-mask FILE haplotype indices (0-based) allowed to contribute strings [all]\n");
		fprintf(fp_help, "    --mate-resolve     an end whose mate has mapq>=--mate-mapq keeps only candidates\n");
		fprintf(fp_help, "                       compatible with that mate, nearest first, BEFORE mapq; no 0.2/0.8 blend\n");
		fprintf(fp_help, "    --mate-mapq INT    standalone mapq that makes a mate unambiguous [%d]\n", opt.mate_mapq);
		fprintf(fp_help, "    --mate-frag INT    typical fragment length: among mate-compatible loci take the one nearest to it,\n");
		fprintf(fp_help, "                       fragment measured with clips and indels; 0 = the shortest span [%d]\n", opt.mate_frag);
		fprintf(fp_help, "    --pair-mapq STR    off | vg: pair score dp0+dp1-(d-mate_frag)^2/(2 mate_sd^2 lambda), mapq = softmax over\n");
		fprintf(fp_help, "                       pairs, ONE mapq min(raw,120)/2 for both ends (as vg giraffe) | vgup: vg, ends only raised [off]\n");
		fprintf(fp_help, "    --mate-sd INT      s.d. of the fragment length for the --pair-mapq prior; 0 = no prior [%d]\n", opt.mate_sd);
		fprintf(fp_help, "    --mate-gap INT     --mate-resolve: an end is unambiguous when its best DP score beats every other\n");
		fprintf(fp_help, "                       candidate by >= INT (a+b = one substitution); 0 = mapq >= --mate-mapq and n_sub 0 [%d]\n", opt.mate_gap);
		fprintf(fp_help, "    --mapq-quotient    collapse candidate alignments by projection onto the\n");
		fprintf(fp_help, "                       base reference instead of by sequence name\n");
		fprintf(fp_help, "    --mapq-prefer-primary  on score ties keep the representative on the base\n");
		fprintf(fp_help, "                       reference, so fewer reads need lift-over\n");
        fprintf(fp_help, "    --parse-haplotype parse haplotypes from VCF to generate real SNP combinations, otherwise use all.\n");
		fprintf(fp_help, "  Mapping:\n");
		fprintf(fp_help, "    -f FLOAT     filter out top FLOAT fraction of repetitive minimizers [%g]\n", opt.mid_occ_frac);
		fprintf(fp_help, "    -g NUM       stop chain enlongation if there are no minimizers in INT-bp [%d]\n", opt.max_gap);
		fprintf(fp_help, "    -G NUM       max intron length (effective with -xsplice; changing -r) [200k]\n");
		fprintf(fp_help, "    -F NUM       max fragment length (effective with -xsr or in the fragment mode) [800]\n");
		fprintf(fp_help, "    -r NUM[,NUM] chaining/alignment bandwidth and long-join bandwidth [%d,%d]\n", opt.bw, opt.bw_long);
		fprintf(fp_help, "    -n INT       minimal number of minimizers on a chain [%d]\n", opt.min_cnt);
		fprintf(fp_help, "    -m INT       minimal chaining score (matching bases minus log gap penalty) [%d]\n", opt.min_chain_score);
//		fprintf(fp_help, "    -T INT       SDUST threshold; 0 to disable SDUST [%d]\n", opt.sdust_thres); // TODO: this option is never used; might be buggy
		fprintf(fp_help, "    -X           skip self and dual mappings (for the all-vs-all mode)\n");
		fprintf(fp_help, "    -p FLOAT     min secondary-to-primary score ratio [%g]\n", opt.pri_ratio);
		fprintf(fp_help, "    -N INT       retain at most INT secondary alignments [%d]\n", opt.best_n);
		fprintf(fp_help, "  Alignment:\n");
		fprintf(fp_help, "    -A INT       matching score [%d]\n", opt.a);
		fprintf(fp_help, "    -B INT       mismatch penalty (larger value for lower divergence) [%d]\n", opt.b);
		fprintf(fp_help, "    -O INT[,INT] gap open penalty [%d,%d]\n", opt.q, opt.q2);
		fprintf(fp_help, "    -E INT[,INT] gap extension penalty; a k-long gap costs min{O1+k*E1,O2+k*E2} [%d,%d]\n", opt.e, opt.e2);
		fprintf(fp_help, "    -z INT[,INT] Z-drop score and inversion Z-drop score [%d,%d]\n", opt.zdrop, opt.zdrop_inv);
		fprintf(fp_help, "    -s INT       minimal peak DP alignment score [%d]\n", opt.min_dp_max);
		fprintf(fp_help, "    -u CHAR      how to find GT-AG. f:transcript strand, b:both strands, n:don't match GT-AG [n]\n");
		fprintf(fp_help, "    -J INT       splice mode. 0: original minimap2 model; 1: miniprot model [1]\n");
		fprintf(fp_help, "    -j FILE      junctions in BED12 to extend *short* RNA-seq alignment []\n");
		fprintf(fp_help, "  Input/Output:\n");
		fprintf(fp_help, "    -a           output in the SAM format (PAF by default)\n");
		fprintf(fp_help, "    -o FILE      output alignments to FILE [stdout]\n");
		fprintf(fp_help, "    -L           write CIGAR with >65535 ops at the CG tag\n");
		fprintf(fp_help, "    -R STR       SAM read group line in a format like '@RG\\tID:foo\\tSM:bar' []\n");
		fprintf(fp_help, "    -c           output CIGAR in PAF\n");
		fprintf(fp_help, "    --cs[=STR]   output the cs tag; STR is 'short' (if absent) or 'long' [none]\n");
		fprintf(fp_help, "    --ds         output the ds tag, which is an extension to cs\n");
		fprintf(fp_help, "    --MD         output the MD tag\n");
		fprintf(fp_help, "    --eqx        write =/X CIGAR operators\n");
		fprintf(fp_help, "    -Y           use soft clipping for supplementary alignments\n");
		fprintf(fp_help, "    -y           copy FASTA/Q comments to output SAM\n");
		fprintf(fp_help, "    -t INT       number of threads [%d]\n", n_threads);
		fprintf(fp_help, "    -K NUM       minibatch size for mapping [500M]\n");
//		fprintf(fp_help, "    -v INT       verbose level [%d]\n", mm_verbose);
		fprintf(fp_help, "    --version    show version number\n");
		fprintf(fp_help, "  Preset:\n");
		fprintf(fp_help, "    -x STR       preset (always applied before other options; see minimap2.1 for details) []\n");
		fprintf(fp_help, "                 - lr:hq - accurate long reads (error rate <1%%) against a reference genome\n");
		fprintf(fp_help, "                 - splice/splice:hq - spliced alignment for long reads/accurate long reads\n");
		fprintf(fp_help, "                 - splice:sr - spliced alignment for short RNA-seq reads\n");
		fprintf(fp_help, "                 - asm5/asm10/asm20 - asm-to-ref mapping, for ~0.1/1/5%% sequence divergence\n");
		fprintf(fp_help, "                 - sr - short reads against a reference\n");
		fprintf(fp_help, "                 - map-pb/map-hifi/map-ont/map-iclr - CLR/HiFi/Nanopore/ICLR vs reference mapping\n");
		fprintf(fp_help, "                 - ava-pb/ava-ont - PacBio CLR/Nanopore read overlap\n");
		fprintf(fp_help, "\nSee `man ./minimap2.1' for detailed description of these and other advanced command-line options.\n");
		return fp_help == stdout? 0 : 1;
	}

	if ((opt.flag & MM_F_SR) && argc - o.ind > 3) {
		fprintf(stderr, "[ERROR] incorrect input: in the sr mode, please specify no more than two query files.\n");
		return 1;
	}
	if (fn_hap_dp && fn_var_alt) { /* both re-score the same mismatches */
		fprintf(stderr, "[ERROR] --hap-dp and --var-alt are mutually exclusive\n");
		return 1;
	}
	if (fn_hap_dp && hap_dp_max < 1) {
		fprintf(stderr, "[ERROR] --hap-dp-max must be >= 1 (the base string counts)\n");
		return 1;
	}
	if (fn_hap_dp_mask && !fn_hap_dp) {
		fprintf(stderr, "[ERROR] --hap-dp-mask needs --hap-dp\n");
		return 1;
	}
	if (rmin_init() != 0) return 1;
	if (rmin_active() && (vcf_with_variants == 0 || fnw == 0)) {
		fprintf(stderr, "[ERROR] --read-min-mode needs --vcf-file-with-variants and -d (index build)\n");
		return 1;
	}
	idx_rdr = mm_idx_reader_open(argv[o.ind], &ipt, fnw);
	if (idx_rdr == 0) {
		fprintf(stderr, "[ERROR] failed to open file '%s': %s\n", argv[o.ind], strerror(errno));
		return 1;
	}
	if (!idx_rdr->is_idx && fnw == 0 && argc - o.ind < 2) {
		fprintf(stderr, "[ERROR] missing input: please specify a query file to map or option -d to keep the index\n");
		mm_idx_reader_close(idx_rdr);
		return 1;
	}
	if (opt.best_n == 0 && (opt.flag&MM_F_CIGAR) && mm_verbose >= 2)
		fprintf(stderr, "[WARNING]\033[1;31m `-N 0' reduces alignment accuracy. Please use --secondary=no to suppress secondary alignments.\033[0m\n");
	while ((mi = mm_idx_reader_read(idx_rdr, n_threads, vcf_with_variants)) != 0) {
		int ret;
		if (opt.flag & MM_F_LIFT) {
			int64_t n_unres = 0;
			int n_lift = mm_lift_build(mi, fn_lift_map, lift_jct, &n_unres);
			if (n_lift < 0) {
				fprintf(stderr, "[ERROR] --lift: could not build the projection table\n");
				mm_idx_destroy(mi); mm_idx_reader_close(idx_rdr); return 1;
			}
			if (mm_verbose >= 3)
				fprintf(stderr, "[M::%s] --lift: %d of %u sequences projectable, %ld unresolved, junction from the sidecar for %ld\n",
						__func__, n_lift, mi->n_seq, (long)n_unres, (long)(mi->lift? mi->lift->n_jct_from_map : 0));
			if (mm_verbose >= 3 && mi->lift && mi->lift->n_multi > 0)   /* only for a lift map with a fifth column */
				fprintf(stderr, "[M::%s] --lift: %ld haplotype contigs with several gaps (%ld further gaps in all)\n",
						__func__, (long)mi->lift->n_multi, (long)mi->lift->n_xg);
		}
		if (fn_var_alt) {
			double t0 = realtime();
			mm_varalt_t *va;
			if (mi->flag & MM_I_NO_SEQ) {
				fprintf(stderr, "[ERROR] --var-alt needs the index with sequences\n");
				mm_idx_destroy(mi); mm_idx_reader_close(idx_rdr); return 1;
			}
			va = mm_varalt_load(mi, fn_var_alt);
			if (va == NULL) {
				fprintf(stderr, "[ERROR] --var-alt: failed to read '%s'\n", fn_var_alt);
				mm_idx_destroy(mi); mm_idx_reader_close(idx_rdr); return 1;
			}
			mi->varalt = va;
			if (mm_verbose >= 3)
				fprintf(stderr, "[M::%s] --var-alt: %ld lines, %ld positions, %ld ALT alleles kept (%.1f MB); unknown sequence %ld, malformed %ld, "
						"ALT = base %ld, beyond end %ld, merged duplicates %ld; loaded in %.2f s\n", __func__,
						(long)va->n_lines, (long)va->n_pos, (long)va->n_alleles, va->n_pos * 5.0 / 1048576.0, (long)va->n_unknown_seq,
						(long)va->n_bad, (long)va->n_eq_ref, (long)va->n_out_of_range, (long)va->n_dup, realtime() - t0);
		}
		if (fn_hap_dp) {
			double t0 = realtime();
			mm_hapdp_t *hd;
			if (mi->flag & MM_I_NO_SEQ) {
				fprintf(stderr, "[ERROR] --hap-dp needs the index with sequences\n");
				mm_idx_destroy(mi); mm_idx_reader_close(idx_rdr); return 1;
			}
			hd = mm_hapdp_load(mi, fn_hap_dp);
			if (hd == NULL) {
				fprintf(stderr, "[ERROR] --hap-dp: failed to read '%s'\n", fn_hap_dp);
				mm_idx_destroy(mi); mm_idx_reader_close(idx_rdr); return 1;
			}
			hd->max_str = hap_dp_max;
			if (fn_hap_dp_mask && mm_hapdp_load_mask(hd, fn_hap_dp_mask) != 0) {
				fprintf(stderr, "[ERROR] --hap-dp-mask: failed to read '%s' (or a haplotype index out of range)\n", fn_hap_dp_mask);
				mm_hapdp_destroy(hd); mm_idx_destroy(mi); mm_idx_reader_close(idx_rdr); return 1;
			}
			mi->hapdp = hd;
			if (mm_verbose >= 3)
				fprintf(stderr, "[M::%s] --hap-dp: %ld lines, %ld SNV records kept, %d haplotypes, max strings %d, mask %s; unknown sequence %ld, "
						"malformed %ld, ALT = base %ld, beyond end %ld, no carrier %ld, merged duplicates %ld, carrier conflicts %ld; loaded in %.2f s\n",
						__func__, (long)hd->n_lines, (long)hd->n_rec, hd->n_hap, hd->max_str, hd->has_allow ? "yes" : "no",
						(long)hd->n_unknown_seq, (long)hd->n_bad, (long)hd->n_eq_ref, (long)hd->n_out_of_range, (long)hd->n_nocarrier,
						(long)hd->n_dup, (long)hd->n_conflict, realtime() - t0);
		}
		if ((opt.flag & MM_F_CIGAR) && (mi->flag & MM_I_NO_SEQ)) {
			fprintf(stderr, "[ERROR] the prebuilt index doesn't contain sequences.\n");
			mm_idx_destroy(mi);
			mm_idx_reader_close(idx_rdr);
			return 1;
		}
		if ((opt.flag & MM_F_OUT_SAM) && idx_rdr->n_parts == 1) {
			if (mm_idx_reader_eof(idx_rdr)) {
				if (opt.split_prefix == 0)
					ret = mm_write_sam_hdr(mi, rg, MM_VERSION, argc, argv);
				else
					ret = mm_write_sam_hdr(0, rg, MM_VERSION, argc, argv);
			} else {
				ret = mm_write_sam_hdr(0, rg, MM_VERSION, argc, argv);
				if (opt.split_prefix == 0 && mm_verbose >= 2)
					fprintf(stderr, "[WARNING]\033[1;31m For a multi-part index, no @SQ lines will be outputted. Please use --split-prefix.\033[0m\n");
			}
			if (ret != 0) {
				mm_idx_destroy(mi);
				mm_idx_reader_close(idx_rdr);
				return 1;
			}
		}
		if (mm_verbose >= 3)
			fprintf(stderr, "[M::%s::%.3f*%.2f] loaded/built the index for %d target sequence(s)\n",
					__func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), mi->n_seq);
		if (argc != o.ind + 1) mm_mapopt_update(&opt, mi);
		if (mm_verbose >= 3) mm_idx_stat(mi);
		rmin_report(stderr);
		if (fn_bed_junc) {
			mm_idx_bed_read(mi, fn_bed_junc, 1);
			if (mi->I == 0 && mm_verbose >= 2)
				fprintf(stderr, "[WARNING] failed to load the junction BED file\n");
		}
		if (fn_bed_jump) {
			mm_idx_jjump_read(mi, fn_bed_jump, MM_JUNC_ANNO, -1);
			if (mi->J == 0 && mm_verbose >= 2)
				fprintf(stderr, "[WARNING] failed to load the jump BED file\n");
		}
		if (fn_bed_pass1) {
			mm_idx_jjump_read(mi, fn_bed_pass1, MM_JUNC_MISC, 5);
			if (mi->J == 0 && mm_verbose >= 2)
				fprintf(stderr, "[WARNING] failed to load the pass-1 jump BED file\n");
		}
		if (fn_spsc) {
			mm_idx_spsc_read2(mi, fn_spsc, mm_max_spsc_bonus(&opt), spsc_scale);
			if (mi->spsc == 0 && mm_verbose >= 2)
				fprintf(stderr, "[WARNING] failed to load the splice score file\n");
		}
		if (alt_list) mm_idx_alt_read(mi, alt_list);
		if (argc - (o.ind + 1) == 0) {
			mm_idx_destroy(mi);
			continue; // no query files
		}
		ret = 0;
		if (!(opt.flag & MM_F_FRAG_MODE)) {
			for (i = o.ind + 1; i < argc; ++i) {
				ret = mm_map_file(mi, argv[i], &opt, n_threads);
				if (ret < 0) break;
			}
		} else {
			ret = mm_map_file_frag(mi, argc - (o.ind + 1), (const char**)&argv[o.ind + 1], &opt, n_threads);
		}
		mm_idx_destroy(mi);
		if (ret < 0) {
			fprintf(stderr, "ERROR: failed to map the query file\n");
			exit(EXIT_FAILURE);
		}
	}
	if (idx_rdr->error || idx_rdr->n_parts == 0) {
		if (idx_rdr->error)
			fprintf(stderr, "[ERROR] could not %s '%s' at part %d: %s\n",
					idx_rdr->is_idx ? "read binary index" : "build reference index",
					argv[o.ind], idx_rdr->n_parts + 1, strerror(idx_rdr->error));
		else fprintf(stderr, "[ERROR] no index or reference sequences read from '%s'\n", argv[o.ind]);
		mm_idx_reader_close(idx_rdr);
		return 1;
	}
	n_parts = idx_rdr->n_parts;
	mm_idx_reader_close(idx_rdr);

	if (opt.split_prefix)
		mm_split_merge(argc - (o.ind + 1), (const char**)&argv[o.ind + 1], &opt, n_parts);

	if (fflush(stdout) == EOF) {
		perror("[ERROR] failed to write the results");
		exit(EXIT_FAILURE);
	}

	if (mm_verbose >= 3) {
		fprintf(stderr, "[M::%s] Version: %s\n", __func__, MM_VERSION);
		fprintf(stderr, "[M::%s] CMD:", __func__);
		for (i = 0; i < argc; ++i)
			fprintf(stderr, " %s", argv[i]);
		fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, realtime() - mm_realtime0, cputime(), peakrss() / 1024.0 / 1024.0 / 1024.0);
	}
	return 0;
}
