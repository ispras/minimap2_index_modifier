/* Unit tests for mm_lift.c.
 *
 * Build:  cc -O1 -Wall -I.. -o test_mm_lift test_mm_lift.c ../mm_lift.c
 * Run:    ./test_mm_lift            (expects chain_sample.txt next to it)
 *
 * The fixtures are real: test/chain_sample.txt is the first 400 records of
 * gen2/fixed2.chain on 172.16.13.35, and test/fai_sample.tsv the matching
 * lines of gen2/fixed2.fasta.fai.  Test 6 re-derives the projection from the
 * contig name + length alone and checks it against the chain's own header
 * coordinates for all 400, which is the claim DESIGN.md rests on.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mm_lift.h"

#define OP_M 0
#define OP_I 1
#define OP_D 2
#define OP_S 4
#define CIG_MK(l, o) ((uint32_t)(l) << 4 | (o))
#define CIG_LEN(c) ((c) >> 4)
#define CIG_OP(c)  ((c) & 0xf)

static int fails = 0, checks = 0;
static const char *OPS = "MIDNSHP=X";

static void cig_str(const uint32_t *c, int32_t n, char *out)
{
	int i, k = 0;
	if (n == 0) { strcpy(out, "<empty>"); return; }
	for (i = 0; i < n; ++i)
		k += sprintf(out + k, "%u%c", (unsigned)CIG_LEN(c[i]), OPS[CIG_OP(c[i])]);
}

static void expect_cigar(const char *what, const uint32_t *got, int32_t n_got,
                         const char *want)
{
	char buf[512];
	++checks;
	cig_str(got, n_got, buf);
	if (strcmp(buf, want) != 0) {
		printf("  FAIL %-46s cigar %s, want %s\n", what, buf, want);
		++fails;
	}
}

static void expect_i64(const char *what, int64_t got, int64_t want)
{
	++checks;
	if (got != want) {
		printf("  FAIL %-46s got %lld, want %lld\n", what,
		       (long long)got, (long long)want);
		++fails;
	}
}

static void expect_int(const char *what, int got, int want)
{
	++checks;
	if (got != want) {
		printf("  FAIL %-46s got %d, want %d\n", what, got, want);
		++fails;
	}
}

/* ---- the two real contigs used throughout ----------------------------------
 * chr1_contig_30884  len 205  chain "101 4<TAB>0"  insertion, ins_len 4
 * chr1_contig_54750  len 201  chain "101 0<TAB>4"  deletion,  del_len 4
 * both with qStart == pos - 101, which is off_a.
 */
static const mm_lift_rec_t INS = { 0, 101, 4, 0, 30783 };  /* chr1_contig_30884 */
static const mm_lift_rec_t DEL = { 0, 101, 0, 4, 54649 };  /* chr1_contig_54750 */

static void t1_parse_name(void)
{
	int32_t cl; int64_t pos;
	puts("[1] name parsing");

	expect_int("chr1_contig_30884 parses",
	           mm_lift_parse_name("chr1_contig_30884", &cl, &pos), 1);
	expect_int("  chrom_len", cl, 4);
	expect_i64("  pos", pos, 30884);

	/* the two-digit chromosome that gen2/contig_generation.py truncated to
	 * "chr1" via contig_header[1:5]; must survive intact here */
	expect_int("chr10_contig_5 parses",
	           mm_lift_parse_name("chr10_contig_5", &cl, &pos), 1);
	expect_int("  chrom_len", cl, 5);
	expect_i64("  pos", pos, 5);

	expect_int("chrX_contig_999 parses",
	           mm_lift_parse_name("chrX_contig_999", &cl, &pos), 1);
	expect_int("  chrom_len", cl, 4);

	expect_int("plain chromosome rejected",
	           mm_lift_parse_name("chr1", &cl, &pos), 0);
	expect_int("no digits rejected",
	           mm_lift_parse_name("chr1_contig_", &cl, &pos), 0);
	expect_int("trailing junk rejected",
	           mm_lift_parse_name("chr1_contig_12x", &cl, &pos), 0);
	expect_int("empty chrom rejected",
	           mm_lift_parse_name("_contig_12", &cl, &pos), 0);
	/* last occurrence wins, so a chromosome whose name embeds the token is
	 * split at the right place */
	expect_int("embedded token: last wins",
	           mm_lift_parse_name("a_contig_b_contig_7", &cl, &pos), 1);
	expect_int("  chrom_len", cl, 10);
	expect_i64("  pos", pos, 7);
}

static void t2_insertion(void)
{
	uint32_t in[8], out[16];
	int32_t n_out, qso, qeo;
	int64_t rso;
	int rc;
	puts("[2] insertion contig chr1_contig_30884 (ins_len 4, off_a 30783)");

	/* entirely inside block A: pure shift, no surgery */
	in[0] = CIG_MK(50, OP_M);
	rc = mm_lift_cigar(&INS, 205, 0, in, 1, 0, 50, 0, &rso, &qso, &qeo, out, 16, &n_out);
	expect_int("block A rc", rc, 0);
	expect_cigar("block A cigar", out, n_out, "50M");
	expect_i64("block A rs", rso, 30783);

	/* entirely inside block B: shift by off_b = off_a - ins_len = 30779 */
	in[0] = CIG_MK(100, OP_M);
	rc = mm_lift_cigar(&INS, 205, 105, in, 1, 0, 100, 0, &rso, &qso, &qeo, out, 16, &n_out);
	expect_int("block B rc", rc, 0);
	expect_cigar("block B cigar", out, n_out, "100M");
	expect_i64("block B rs", rso, 30884);   /* off_b 30779 + 105 */

	/* THE WORKED EXAMPLE: 150 bp read spanning flank / ALT / flank.
	 * contig [51, 201): 50 bases of block A, the 4 novel ALT bases, 96 of
	 * block B.  The 4 novel bases have no base-reference counterpart, so they
	 * must become an insertion -- NOT a deletion, and not silently dropped. */
	in[0] = CIG_MK(150, OP_M);
	rc = mm_lift_cigar(&INS, 205, 51, in, 1, 0, 150, 0, &rso, &qso, &qeo, out, 16, &n_out);
	expect_int("junction rc", rc, 0);
	expect_cigar("junction cigar", out, n_out, "50M4I96M");
	expect_i64("junction rs", rso, 30834);        /* off_a 30783 + 51 */
	expect_int("junction qs unchanged", qso, 0);
	expect_int("junction qe unchanged", qeo, 150);
	/* G2: query length preserved */
	expect_int("junction query len", 50 + 4 + 96, qeo - qso);

	/* read that stops exactly at the junction: no insertion, no deletion */
	in[0] = CIG_MK(50, OP_M);
	rc = mm_lift_cigar(&INS, 205, 51, in, 1, 0, 50, 0, &rso, &qso, &qeo, out, 16, &n_out);
	expect_int("stops at junction rc", rc, 0);
	expect_cigar("stops at junction cigar", out, n_out, "50M");

	/* read starting inside the novel ALT bases: those query bases cannot be
	 * placed, so they become clipping and qs moves */
	in[0] = CIG_MK(98, OP_M);           /* contig [103, 201) = 2 novel + 96 B */
	rc = mm_lift_cigar(&INS, 205, 103, in, 1, 0, 98, 0, &rso, &qso, &qeo, out, 16, &n_out);
	expect_int("starts in gap rc", rc, 0);
	expect_cigar("starts in gap cigar", out, n_out, "96M");
	expect_i64("starts in gap rs", rso, 30884);
	expect_int("starts in gap qs moved", qso, 2);
	expect_int("starts in gap qe", qeo, 98);
}

static void t3_deletion(void)
{
	uint32_t in[8], out[16];
	int32_t n_out, qso, qeo;
	int64_t rso;
	int rc;
	puts("[3] deletion contig chr1_contig_54750 (del_len 4, off_a 54649)");

	/* entirely inside block A: no deletion invented */
	in[0] = CIG_MK(50, OP_M);
	rc = mm_lift_cigar(&DEL, 201, 0, in, 1, 0, 50, 0, &rso, &qso, &qeo, out, 16, &n_out);
	expect_int("block A rc", rc, 0);
	expect_cigar("block A cigar (no D)", out, n_out, "50M");
	expect_i64("block A rs", rso, 54649);

	/* entirely inside block B: shifted, still no deletion */
	in[0] = CIG_MK(100, OP_M);
	rc = mm_lift_cigar(&DEL, 201, 101, in, 1, 0, 100, 0, &rso, &qso, &qeo, out, 16, &n_out);
	expect_int("block B rc", rc, 0);
	expect_cigar("block B cigar (no D)", out, n_out, "100M");
	expect_i64("block B rs", rso, 54754);   /* off_b 54653 + 101 */

	/* THE WORKED EXAMPLE: 150 bp read spanning the deletion junction.
	 * In contig space it is a clean 150M -- the contig omits the 4 deleted
	 * bases.  In base-reference space those 4 bases ARE skipped, so exactly
	 * one 4D belongs at the junction.  This is the CIGAR levioSAM2 also
	 * produces here, and it is correct; the measured pathology is the same
	 * op emitted when the contig sequence did NOT omit the bases. */
	in[0] = CIG_MK(150, OP_M);
	rc = mm_lift_cigar(&DEL, 201, 51, in, 1, 0, 150, 0, &rso, &qso, &qeo, out, 16, &n_out);
	expect_int("junction rc", rc, 0);
	expect_cigar("junction cigar", out, n_out, "50M4D100M");
	expect_i64("junction rs", rso, 54700);  /* off_a 54649 + 51 */
	expect_int("junction query len", 150, qeo - qso);

	/* G5, the anti-pathology test: a read that does NOT cross the junction
	 * must never receive the deletion, however close it sits. */
	in[0] = CIG_MK(51, OP_M);           /* contig [50, 101), ends AT the junction */
	rc = mm_lift_cigar(&DEL, 201, 50, in, 1, 0, 51, 0, &rso, &qso, &qeo, out, 16, &n_out);
	expect_int("abuts junction rc", rc, 0);
	expect_cigar("abuts junction: no D", out, n_out, "51M");
	expect_i64("abuts junction rs", rso, 54699);

	in[0] = CIG_MK(50, OP_M);           /* contig [101, 151), starts AT it */
	rc = mm_lift_cigar(&DEL, 201, 101, in, 1, 0, 50, 0, &rso, &qso, &qeo, out, 16, &n_out);
	expect_int("starts at junction rc", rc, 0);
	expect_cigar("starts at junction: no D", out, n_out, "50M");
	expect_i64("starts at junction rs", rso, 54754);

	/* the deletion op appears exactly once even for a long span */
	in[0] = CIG_MK(200, OP_M);
	rc = mm_lift_cigar(&DEL, 201, 0, in, 1, 0, 200, 0, &rso, &qso, &qeo, out, 16, &n_out);
	expect_int("full contig rc", rc, 0);
	expect_cigar("full contig cigar", out, n_out, "101M4D99M");
}

static void t4_existing_indels(void)
{
	uint32_t in[8], out[16];
	int32_t n_out, qso, qeo;
	int64_t rso;
	int rc;
	puts("[4] alignments that already contain indels");

	/* read with its own 3 bp insertion before a deletion junction */
	in[0] = CIG_MK(30, OP_M); in[1] = CIG_MK(3, OP_I); in[2] = CIG_MK(117, OP_M);
	/* contig [51, 198): 50 in A, 97 in B */
	rc = mm_lift_cigar(&DEL, 201, 51, in, 3, 0, 150, 0, &rso, &qso, &qeo, out, 16, &n_out);
	expect_int("own I then junction rc", rc, 0);
	expect_cigar("own I then junction", out, n_out, "30M3I20M4D97M");
	expect_i64("own I then junction rs", rso, 54700);

	/* read with its own deletion that ends right at the junction: the
	 * junction D must merge with it, not sit beside it (G3) */
	in[0] = CIG_MK(40, OP_M); in[1] = CIG_MK(11, OP_D); in[2] = CIG_MK(110, OP_M);
	/* contig [50, 211) -> A [50,101)=51 ... but 40M takes [50,90), 11D [90,101),
	 * then 110M from 101 in B */
	rc = mm_lift_cigar(&DEL, 201, 50, in, 3, 0, 150, 0, &rso, &qso, &qeo, out, 16, &n_out);
	expect_int("own D abutting junction rc", rc, 0);
	expect_cigar("own D merges with junction D", out, n_out, "40M15D110M");

	/* insertion contig: a read whose own deletion sits inside the novel ALT
	 * bases -- those contig bases do not exist in the base reference either,
	 * so the op must vanish entirely, not become a base-reference deletion */
	in[0] = CIG_MK(50, OP_M); in[1] = CIG_MK(4, OP_D); in[2] = CIG_MK(96, OP_M);
	/* contig [51,101) A, [101,105) the gap, [105,201) B */
	rc = mm_lift_cigar(&INS, 205, 51, in, 3, 0, 146, 0, &rso, &qso, &qeo, out, 16, &n_out);
	expect_int("D inside gap rc", rc, 0);
	expect_cigar("D inside gap vanishes", out, n_out, "146M");
	expect_i64("D inside gap rs", rso, 30834);
}

static void t5_edges(void)
{
	uint32_t in[8], out[16];
	int32_t n_out, qso, qeo;
	int64_t rso;
	int rc;
	puts("[5] error and edge cases");

	/* alignment entirely inside the novel ALT bases -> no representation.
	 * Needs a gap wide enough to contain it, which this reference's own
	 * generator forbids (ALT <= 100 bp, reads 150 bp); constructed here only
	 * to prove the code reports it rather than emitting nonsense. */
	{
		mm_lift_rec_t big = { 0, 101, 120, 0, 1000 };
		in[0] = CIG_MK(50, OP_M);
		rc = mm_lift_cigar(&big, 321, 130, in, 1, 0, 50, 0, &rso, &qso, &qeo, out, 16, &n_out);
		expect_int("wholly in gap -> -2", rc, -2);
	}

	/* capacity too small -> -1 and the required op count */
	in[0] = CIG_MK(150, OP_M);
	rc = mm_lift_cigar(&INS, 205, 51, in, 1, 0, 150, 0, &rso, &qso, &qeo, out, 1, &n_out);
	expect_int("too small -> -1", rc, -1);
	expect_int("required op count", n_out, 3);

	/* identity record (a base sequence) is never passed here, but a contig
	 * with neither insertion nor deletion must be a pure shift */
	{
		mm_lift_rec_t flat = { 0, 101, 0, 0, 500 };
		in[0] = CIG_MK(150, OP_M);
		rc = mm_lift_cigar(&flat, 201, 10, in, 1, 0, 150, 0, &rso, &qso, &qeo, out, 16, &n_out);
		expect_int("flat rc", rc, 0);
		expect_cigar("flat cigar", out, n_out, "150M");
		expect_i64("flat rs", rso, 510);
	}
}

/* ---- test 6: re-derive the map from (name, len) and check it against the
 * real chain's own header and block lines, for all 400 sample records. ---- */
static void t6_against_real_chain(void)
{
	FILE *fp;
	char line[1024];
	long n_rec = 0, n_off_ok = 0, n_geom_ok = 0;
	long n_blockB_bad = 0;
	char name[256]; long long tsize, tstart, tend, qstart, qend, qsize;
	char qname[128];
	int have_hdr = 0;
	long long first = 0, dt = 0, dq = 0;

	puts("[6] re-derived map vs. 400 real records of gen2/fixed2.chain");
	fp = fopen("chain_sample.txt", "r");
	if (fp == NULL) fp = fopen("test/chain_sample.txt", "r");
	if (fp == NULL) { printf("  SKIP (chain_sample.txt not found)\n"); return; }

	while (fgets(line, sizeof line, fp)) {
		if (strncmp(line, "chain ", 6) == 0) {
			long long score;
			if (sscanf(line, "chain %lld %255s %lld %*s %lld %lld %127s %lld %*s %lld %lld",
			           &score, name, &tsize, &tstart, &tend,
			           qname, &qsize, &qstart, &qend) == 9) {
				have_hdr = 1; first = dt = dq = 0;
			}
			continue;
		}
		if (!have_hdr) continue;
		{	/* separate temporaries: a partial match on the closing line would
			 * otherwise clobber `first` before we get to read it */
			long long a, b, c;
			if (sscanf(line, "%lld %lld %lld", &a, &b, &c) == 3) {
				first = a; dt = b; dq = c;
				continue;
			}
		}
		if (line[0] == '\n' || line[0] == '\0') continue;
		{
			/* the single-number line closes the record */
			long long last;
			if (sscanf(line, "%lld", &last) != 1) continue;
			have_hdr = 0;
			if (strstr(name, "_contig_") == NULL) continue;  /* chromosome record */
			++n_rec;

			{
				int32_t cl; int64_t pos;
				mm_lift_rec_t rec;
				if (!mm_lift_parse_name(name, &cl, &pos)) {
					printf("  FAIL name did not parse: %s\n", name); ++fails; continue;
				}
				++checks;
				rec.base_rid = 0;
				rec.jct = 101;
				rec.ins_len = (int32_t)dt;   /* chain's target gap = novel ALT bases */
				rec.del_len = (int32_t)dq;   /* chain's query gap  = deleted REF bases */
				rec.off_a = pos - 101;

				/* claim 1: off_a derived from the NAME alone equals the chain's
				 * own qStart.  This is what lets the index carry no coordinates. */
				if (rec.off_a == qstart) ++n_off_ok;

				/* claim 2: jct is 101 for every record */
				if (first == 101) ++n_geom_ok;

				/* claim 3 (the defect): the chain's closing block length is NOT
				 * the true 100 bp right flank.  Insertions get 99, deletions
				 * get 98+|REF|.  Counting how often it disagrees. */
				if (last != 100) ++n_blockB_bad;
			}
		}
	}
	fclose(fp);

	printf("  records parsed                     : %ld\n", n_rec);
	expect_i64("  off_a from name == chain qStart", n_off_ok, n_rec);
	expect_i64("  junction is 101 in every record", n_geom_ok, n_rec);
	/* Not a failure of our code -- a measurement of the chain's defect.
	 * log/t1pilot/chainfix.log step [3] found 0.0342 % of all 570,065 chain
	 * records internally consistent, i.e. only the 195 chromosome records. */
	printf("  chain closing block != 100 bp      : %ld of %ld  (the Д2 defect;\n"
	       "                                       our map never reads it)\n",
	       n_blockB_bad, n_rec);
}


/* ---- test 7: reverse-strand clip orientation -------------------------------
 * Regression for a defect found by code review, not by any length check.
 * mm_lift_cigar() reports clipping in CIGAR (reference) order, but qs/qe live
 * in forward-read order.  write_sam_cigar() prints clip_len[0] before the
 * CIGAR and computes it as (rev ? qlen-qe : qs) -- format.c:490-491 -- so for
 * a reverse alignment the leading clip must come out of qe, not qs.
 *
 * Applying it to the wrong end still leaves
 *     clip_len[0] + cigar_query + clip_len[1] == qlen,
 * so gate G4 (CIGAR query length equals SEQ length) passes either way.  Only
 * this test distinguishes them.
 * ------------------------------------------------------------------------ */
static void t7_reverse_clip(void)
{
	uint32_t in[4], out[16];
	int32_t n_out, qso_f, qeo_f, qso_r, qeo_r;
	int64_t rso;
	int rc;
	puts("[7] reverse-strand clip lands on the other end");

	/* 98 bp alignment starting at contig offset 103, i.e. 2 bases inside the
	 * 4 novel ALT bases of chr1_contig_30884.  Those 2 query bases cannot be
	 * placed and become clipping. */
	in[0] = CIG_MK(98, OP_M);

	rc = mm_lift_cigar(&INS, 205, 103, in, 1, 0, 98, 0,
	                   &rso, &qso_f, &qeo_f, out, 16, &n_out);
	expect_int("forward rc", rc, 0);
	expect_cigar("forward cigar", out, n_out, "96M");
	expect_int("forward qs advances by 2", qso_f, 2);
	expect_int("forward qe unchanged", qeo_f, 98);

	rc = mm_lift_cigar(&INS, 205, 103, in, 1, 0, 98, 1,
	                   &rso, &qso_r, &qeo_r, out, 16, &n_out);
	expect_int("reverse rc", rc, 0);
	expect_cigar("reverse cigar", out, n_out, "96M");
	expect_int("reverse qs unchanged", qso_r, 0);
	expect_int("reverse qe retreats by 2", qeo_r, 96);

	/* both orientations keep the query length, which is exactly why a length
	 * check cannot tell them apart */
	expect_int("forward query span", qeo_f - qso_f, 96);
	expect_int("reverse query span", qeo_r - qso_r, 96);
	/* ...but they are genuinely different placements */
	++checks;
	if (qso_f == qso_r && qeo_f == qeo_r) {
		printf("  FAIL reverse orientation is not distinguished from forward\n");
		++fails;
	}
}

/* ---- 8. same-locus test: shared diagonal, sign-sensitive ------------ */
static mm_reg1_t mk_reg(int rid, int rev, int rs, int qs, int qe, const uint32_t *cig, int n, mm_extra_t **keep)
{
	mm_reg1_t r; memset(&r, 0, sizeof r);
	r.rid = rid; r.rev = rev; r.rs = rs; r.qs = qs; r.qe = qe;
	if (cig) {
		int i, span = 0;
		mm_extra_t *p = (mm_extra_t*)calloc(1, sizeof(mm_extra_t) + n * 4);
		p->n_cigar = n; memcpy(p->cigar, cig, n * 4);
		for (i = 0; i < n; ++i) if (CIG_OP(cig[i]) == OP_M || CIG_OP(cig[i]) == OP_D) span += CIG_LEN(cig[i]);
		r.p = p; r.re = rs + span; *keep = p;
	} else r.re = rs + (qe - qs);
	return r;
}

static void t8_same_locus(void)
{
	mm_extra_t *pa, *pb, *pc, *pd, *pe, *pf;
	uint32_t c150[] = { CIG_MK(150, OP_M) };
	uint32_t cins[] = { CIG_MK(50, OP_M), CIG_MK(4, OP_I), CIG_MK(96, OP_M) };   /* projected contig read, insertion */
	uint32_t cdel[] = { CIG_MK(50, OP_M), CIG_MK(4, OP_D), CIG_MK(100, OP_M) };  /* projected contig read, deletion */
	uint32_t cclip[] = { CIG_MK(96, OP_M) };                                     /* base placement anchored right, 54 soft-clipped */
	mm_reg1_t a, b, c, d, e, f;
	puts("[8] same locus: shared diagonal on a common query interval");
	/* insertion contig read projected vs the same read placed on the base with the same insertion */
	a = mk_reg(0, 0, 1000, 0, 150, cins, 3, &pa);
	b = mk_reg(0, 0, 1000, 0, 150, cins, 3, &pb);
	expect_int("identical projected alignments", mm_lift_same_locus(&a, &b), 1);
	/* base placement that anchors only the right block: shares the block-B diagonal */
	c = mk_reg(0, 0, 1050, 54, 150, cclip, 1, &pc);       /* diag = 1050-54 = 996 = (1000+50) - 54  */
	expect_int("right-anchored base placement shares block B", mm_lift_same_locus(&a, &c), 1);
	/* neighbouring tandem copy: shifted by the period, no shared diagonal */
	d = mk_reg(0, 0, 1040, 0, 150, c150, 1, &pd);
	expect_int("tandem neighbour +40 is NOT the same locus", mm_lift_same_locus(&a, &d), 0);
	d.rs = 960; d.re = 1110;
	expect_int("tandem neighbour -40 is NOT the same locus", mm_lift_same_locus(&a, &d), 0);
	/* same coordinates, other strand */
	e = mk_reg(0, 1, 1000, 0, 150, cins, 3, &pe);
	expect_int("other strand is NOT the same locus", mm_lift_same_locus(&a, &e), 0);
	{	/* REVERSE strand, asymmetric clip: the projected read (qs 0, qe 151, 5D inside) and the base
		 * placement clipped by 4 at the read's end (qs 0, qe 147, starts 9 bp to the right).  In
		 * reference order the read runs backwards, so the base placement's read offset at its
		 * reference-left end is qlen - qe = 4, not qs = 0: the two share the block-B anti-diagonal.
		 * (The earlier diag_segs used qs for both strands and called these two different loci --
		 * two thirds of the mr:i marks on L1aug.) */
		uint32_t crd[] = { CIG_MK(50, OP_M), CIG_MK(5, OP_D), CIG_MK(101, OP_M) };
		uint32_t crc[] = { CIG_MK(147, OP_M) };
		mm_extra_t *p1, *p2, *p3;
		mm_reg1_t g1 = mk_reg(0, 1, 13547254, 0, 151, crd, 3, &p1);   /* rev, 5D at read pos 101 (in ref order after 50M) */
		mm_reg1_t g2 = mk_reg(0, 1, 13547263, 0, 147, crc, 1, &p2);   /* rev, clipped 4 at qe, same place, no indel */
		mm_reg1_t g3 = mk_reg(0, 1, 13547263, 4, 151, crc, 1, &p3);   /* rev, clipped 4 at qs instead: offset differs -> other diagonal */
		expect_int("reverse strand, clip at the read end: same locus", mm_lift_same_locus(&g1, &g2), 1);
		expect_int("reverse strand, clip at the read start: NOT the same locus", mm_lift_same_locus(&g1, &g3), 0);
		free(p1); free(p2); free(p3);
	}
	/* other sequence */
	f = mk_reg(1, 0, 1000, 0, 150, cins, 3, &pf);
	expect_int("other rid is NOT the same locus", mm_lift_same_locus(&a, &f), 0);
	/* deletion: block B diagonal differs from block A by +4; a base read on block A only */
	{
		mm_extra_t *pg, *ph;
		mm_reg1_t g = mk_reg(0, 0, 1000, 0, 150, cdel, 3, &pg);
		uint32_t c50[] = { CIG_MK(50, OP_M) }, c100[] = { CIG_MK(100, OP_M) };
		mm_reg1_t h = mk_reg(0, 0, 1000, 0, 50, c50, 1, &ph);
		expect_int("deletion contig read vs base read on block A", mm_lift_same_locus(&g, &h), 1);
		h.rs = 1004; h.re = 1054;                          /* diag 1004 on q in [0,50): block B's diagonal, block A's query range */
		expect_int("block-B diagonal but query interval disjoint", mm_lift_same_locus(&g, &h), 0);
		free(ph);
		h = mk_reg(0, 0, 1054, 50, 150, c100, 1, &ph);     /* diag 1004 on q in [50,150): block B */
		expect_int("base read sharing block B of a deletion", mm_lift_same_locus(&g, &h), 1);
		free(pg); free(ph);
	}
	/* no CIGAR: end points only */
	{
		mm_reg1_t x, y; mm_extra_t *dummy = 0;
		x = mk_reg(0, 0, 500, 0, 100, 0, 0, &dummy);
		y = mk_reg(0, 0, 520, 20, 100, 0, 0, &dummy);
		expect_int("no-CIGAR regs on one diagonal", mm_lift_same_locus(&x, &y), 1);
		y.rs = 521; y.re = 601;
		expect_int("no-CIGAR regs off by one", mm_lift_same_locus(&x, &y), 0);
	}
	free(pa); free(pb); free(pc); free(pd); free(pe); free(pf);
}

int main(void)
{
	puts("mm_lift unit tests");
	puts("==================");
	t1_parse_name();
	t2_insertion();
	t3_deletion();
	t4_existing_indels();
	t5_edges();
	t6_against_real_chain();
	t7_reverse_clip();
	t8_same_locus();
	printf("\n%d checks, %d failures\n", checks, fails);
	return fails ? 1 : 0;
}
