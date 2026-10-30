/* Unit tests for mm_lift.c.
 *
 * Build:  cc -O1 -Wall -I.. -o test_mm_lift test_mm_lift.c ../mm_lift.c
 * Run:    ./test_mm_lift            (expects chain_sample.txt next to it)
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

	/* a two-digit chromosome name must survive intact */
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
	/* query length preserved */
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
	 * 150M on the contig, which omits the 4 deleted bases, so exactly one
	 * 4D belongs at the junction on the base reference. */
	in[0] = CIG_MK(150, OP_M);
	rc = mm_lift_cigar(&DEL, 201, 51, in, 1, 0, 150, 0, &rso, &qso, &qeo, out, 16, &n_out);
	expect_int("junction rc", rc, 0);
	expect_cigar("junction cigar", out, n_out, "50M4D100M");
	expect_i64("junction rs", rso, 54700);  /* off_a 54649 + 51 */
	expect_int("junction query len", 150, qeo - qso);

	/* a read that does not cross the junction never gets the deletion */
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

	/* an alignment that starts in block B with several reference-consuming
	 * ops (split by the read's own insertion) gets no deletion */
	in[0] = CIG_MK(20, OP_M); in[1] = CIG_MK(3, OP_I); in[2] = CIG_MK(30, OP_M);
	rc = mm_lift_cigar(&DEL, 201, 110, in, 3, 0, 53, 0, &rso, &qso, &qeo, out, 16, &n_out);
	expect_int("starts in B, own I rc", rc, 0);
	expect_cigar("starts in B, own I: no phantom D", out, n_out, "20M3I30M");
	expect_i64("starts in B, own I rs", rso, 54763);   /* off_b 54653 + 110 */
	in[0] = CIG_MK(20, OP_M); in[1] = CIG_MK(2, OP_D); in[2] = CIG_MK(30, OP_M);
	rc = mm_lift_cigar(&DEL, 201, 110, in, 3, 0, 50, 0, &rso, &qso, &qeo, out, 16, &n_out);
	expect_cigar("starts in B, own D: stays 2D", out, n_out, "20M2D30M");

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
	 * junction D must merge with it */
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

	/* alignment entirely inside the novel ALT bases -> no representation */
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

				/* off_a derived from the name alone equals the chain's qStart */
				if (rec.off_a == qstart) ++n_off_ok;

				/* jct is 101 for every record */
				if (first == 101) ++n_geom_ok;

				/* the chain's closing block length is not the 100 bp right flank;
				 * count how often it disagrees */
				if (last != 100) ++n_blockB_bad;
			}
		}
	}
	fclose(fp);

	printf("  records parsed                     : %ld\n", n_rec);
	expect_i64("  off_a from name == chain qStart", n_off_ok, n_rec);
	expect_i64("  junction is 101 in every record", n_geom_ok, n_rec);
	/* informational: a property of the chain file, not a failure */
	printf("  chain closing block != 100 bp      : %ld of %ld  (a known chain defect;\n"
	       "                                       our map never reads it)\n",
	       n_blockB_bad, n_rec);
}


/* ---- test 7: reverse-strand clip orientation -------------------------------
 * mm_lift_cigar() reports clipping in CIGAR (reference) order, but qs/qe are
 * in forward-read order, so for a reverse alignment the leading clip must
 * come out of qe, not qs.  A query-length check cannot tell the two apart.
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
		 * reference-left end is qlen - qe = 4, not qs = 0: the two share the block-B anti-diagonal. */
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

/* ===========================================================================
 * MULTI-GAP (haplotype) contigs.
 *
 * H2: insertion A then deletion B on one haplotype, 40 bp apart:
 *   p1 = 1101, left flank w = 100 -> off_a = p1 - 101 = 1000
 *   A: REF 1 base, ALT 5 bases  -> gap 0: jct 101, ins 4, del 0
 *   B: p2 = 1141, REF 7, ALT 1  -> anchor at contig 100+5+(1141-1101-1) = 144,
 *                                  gap 1: jct 145, ins 0, del 6
 *   contig length 245 = 100 + 5 + 39 + 1 + 100
 * Block offsets: block0 1000, block1 996, block2 1002 (0-based base).
 * ======================================================================== */
static const mm_lift_gap_t H2_X[] = { { 145, 0, 6 } };
static const mm_lift_rec_t H2 = { 0, 101, 4, 0, 1000, 1, H2_X };
static const mm_lift_rec_t H2_A = { 0, 101, 4, 0, 1000, 0, 0 };   /* the single contig of A, same p1 */

static void t9_parse_hap(void)
{
	int32_t cl; int64_t pos;
	puts("[9] haplotype contig names");
	expect_int("chr21_contig_123_h1 parses", mm_lift_parse_name("chr21_contig_123_h1", &cl, &pos), 1);
	expect_int("  chrom_len", cl, 5);
	expect_i64("  pos is the number before _h", pos, 123);
	expect_int("chr21_contig_7_h16 parses", mm_lift_parse_name("chr21_contig_7_h16", &cl, &pos), 1);
	expect_i64("  pos", pos, 7);
	expect_int("_h without digits rejected", mm_lift_parse_name("chr21_contig_123_h", &cl, &pos), 0);
	expect_int("_hx rejected", mm_lift_parse_name("chr21_contig_123_hx", &cl, &pos), 0);
	expect_int("_h1x rejected", mm_lift_parse_name("chr21_contig_123_h1x", &cl, &pos), 0);
	expect_int("_x1 rejected", mm_lift_parse_name("chr21_contig_123_x1", &cl, &pos), 0);
	expect_int("no pos before _h rejected", mm_lift_parse_name("chr21_contig__h1", &cl, &pos), 0);
	expect_int("_h1_h2 rejected", mm_lift_parse_name("chr21_contig_5_h1_h2", &cl, &pos), 0);
}

static void t10_two_gaps(void)
{
	uint32_t in[8], out[32];
	int32_t n_out, qso, qeo;
	int64_t rso;
	int rc;
	puts("[10] two gaps: insertion A then deletion B");

	expect_int("n_gaps", mm_lift_n_gaps(&H2), 2);
	expect_i64("block0 off", mm_lift_block_off(&H2, 0), 1000);
	expect_i64("block1 off == off_b", mm_lift_block_off(&H2, 1), mm_lift_off_b(&H2));
	expect_i64("block2 off", mm_lift_block_off(&H2, 2), 1002);
	expect_int("zone of 100 (A anchor)", mm_lift_zone(&H2, 100), 0);
	expect_int("zone of 101 (A novel)", mm_lift_zone(&H2, 101), 1);
	expect_int("zone of 105 (block 1)", mm_lift_zone(&H2, 105), 2);
	expect_int("zone of 144 (B anchor)", mm_lift_zone(&H2, 144), 2);
	expect_int("zone of 145 (block 2; gap 1 has no novel bases)", mm_lift_zone(&H2, 145), 4);
	expect_i64("pos(100) = p1-1", mm_lift_pos(&H2, 100), 1100);
	expect_i64("pos(102) = junction of A", mm_lift_pos(&H2, 102), 1101);
	expect_i64("pos(105) = p1", mm_lift_pos(&H2, 105), 1101);
	expect_i64("pos(144) = p2-1", mm_lift_pos(&H2, 144), 1140);
	expect_i64("pos(145) = p2+6 (after the 6 deleted)", mm_lift_pos(&H2, 145), 1147);

	/* read across both junctions */
	in[0] = CIG_MK(150, OP_M);
	rc = mm_lift_cigar(&H2, 245, 51, in, 1, 0, 150, 0, &rso, &qso, &qeo, out, 32, &n_out);
	expect_int("across both rc", rc, 0);
	expect_cigar("across both cigar", out, n_out, "50M4I40M6D56M");
	expect_i64("across both rs", rso, 1051);
	expect_int("across both qs", qso, 0);
	expect_int("across both qe", qeo, 150);

	/* only the second junction */
	in[0] = CIG_MK(90, OP_M);
	rc = mm_lift_cigar(&H2, 245, 110, in, 1, 0, 90, 0, &rso, &qso, &qeo, out, 32, &n_out);
	expect_int("second only rc", rc, 0);
	expect_cigar("second only cigar", out, n_out, "35M6D55M");
	expect_i64("second only rs", rso, 1106);

	/* only the first junction, stops before B: no D */
	in[0] = CIG_MK(89, OP_M);
	rc = mm_lift_cigar(&H2, 245, 51, in, 1, 0, 89, 0, &rso, &qso, &qeo, out, 32, &n_out);
	expect_int("first only rc", rc, 0);
	expect_cigar("first only cigar", out, n_out, "50M4I35M");

	/* ends exactly at B's junction (contig 145): still no D */
	in[0] = CIG_MK(94, OP_M);
	rc = mm_lift_cigar(&H2, 245, 51, in, 1, 0, 94, 0, &rso, &qso, &qeo, out, 32, &n_out);
	expect_int("abuts second junction rc", rc, 0);
	expect_cigar("abuts second junction: no D", out, n_out, "50M4I40M");

	/* starts exactly at B's junction: no D, shifted by block 2 */
	in[0] = CIG_MK(80, OP_M);
	rc = mm_lift_cigar(&H2, 245, 145, in, 1, 0, 80, 0, &rso, &qso, &qeo, out, 32, &n_out);
	expect_int("starts at second junction rc", rc, 0);
	expect_cigar("starts at second junction", out, n_out, "80M");
	expect_i64("starts at second junction rs", rso, 1147);

	/* starts inside A's novel bases: 2 clipped, then B's deletion */
	in[0] = CIG_MK(120, OP_M);
	rc = mm_lift_cigar(&H2, 245, 103, in, 1, 0, 120, 0, &rso, &qso, &qeo, out, 32, &n_out);
	expect_int("starts in gap 0 rc", rc, 0);
	expect_cigar("starts in gap 0", out, n_out, "40M6D78M");
	expect_i64("starts in gap 0 rs", rso, 1101);
	expect_int("starts in gap 0 qs", qso, 2);
	expect_int("starts in gap 0 qe", qeo, 120);
	/* the same on the reverse strand: the clip lands at qe */
	rc = mm_lift_cigar(&H2, 245, 103, in, 1, 0, 120, 1, &rso, &qso, &qeo, out, 32, &n_out);
	expect_int("starts in gap 0 rev qs", qso, 0);
	expect_int("starts in gap 0 rev qe", qeo, 118);

	/* ends inside A's novel bases: tail clip, no B */
	in[0] = CIG_MK(43, OP_M);
	rc = mm_lift_cigar(&H2, 245, 60, in, 1, 0, 43, 0, &rso, &qso, &qeo, out, 32, &n_out);
	expect_int("ends in gap 0 rc", rc, 0);
	expect_cigar("ends in gap 0", out, n_out, "41M");
	expect_int("ends in gap 0 qe", qeo, 41);

	/* the read's own insertion right at B's junction, then block 2 */
	in[0] = CIG_MK(94, OP_M); in[1] = CIG_MK(3, OP_I); in[2] = CIG_MK(53, OP_M);
	rc = mm_lift_cigar(&H2, 245, 51, in, 3, 0, 150, 0, &rso, &qso, &qeo, out, 32, &n_out);
	expect_int("own I at second junction rc", rc, 0);
	expect_cigar("own I at second junction", out, n_out, "50M4I40M3I6D53M");

	/* the read's own deletion abutting B's junction merges with B's D */
	in[0] = CIG_MK(80, OP_M); in[1] = CIG_MK(10, OP_D); in[2] = CIG_MK(70, OP_M);
	rc = mm_lift_cigar(&H2, 245, 55, in, 3, 0, 150, 0, &rso, &qso, &qeo, out, 32, &n_out);
	/* contig [55,135) 80M: 46 block0, 4 novel, 30 block1; 10D [135,145); 70M block2 */
	expect_int("own D abutting second junction rc", rc, 0);
	expect_cigar("own D merges with B's D", out, n_out, "46M4I30M16D70M");

	/* capacity: 150M across both gaps needs 5 ops */
	in[0] = CIG_MK(150, OP_M);
	rc = mm_lift_cigar(&H2, 245, 51, in, 1, 0, 150, 0, &rso, &qso, &qeo, out, 2, &n_out);
	expect_int("too small -> -1", rc, -1);
	expect_int("required op count", n_out, 5);
}

/* three gaps, one of them complex (insertion and deletion at once) */
static void t11_three_gaps(void)
{
	/* block0 [0,101) off 5000; gap0 jct 101 ins 3 del 2 (complex); block1 [104,130);
	 * gap1 jct 130 ins 0 del 5; block2 [130,150); gap2 jct 150 ins 7 del 0; block3 [157,257) */
	static const mm_lift_gap_t X[] = { { 130, 0, 5 }, { 150, 7, 0 } };
	static const mm_lift_rec_t R = { 0, 101, 3, 2, 5000, 2, X };
	uint32_t in[4], out[32];
	int32_t n_out, qso, qeo;
	int64_t rso;
	int rc;
	puts("[11] three gaps, one complex");
	expect_int("n_gaps", mm_lift_n_gaps(&R), 3);
	expect_i64("block1 off", mm_lift_block_off(&R, 1), 4999);
	expect_i64("block2 off", mm_lift_block_off(&R, 2), 5004);
	expect_i64("block3 off", mm_lift_block_off(&R, 3), 4997);
	in[0] = CIG_MK(150, OP_M);   /* contig [90, 240) */
	rc = mm_lift_cigar(&R, 257, 90, in, 1, 0, 150, 0, &rso, &qso, &qeo, out, 32, &n_out);
	expect_int("across three rc", rc, 0);
	/* 11 block0, 3 novel (I), D2, 26 block1, D5, 20 block2, 7 novel (I), 83 block3 */
	expect_cigar("across three cigar", out, n_out, "11M3I2D26M5D20M7I83M");
	expect_i64("across three rs", rso, 5090);
	/* ends inside gap 2's novel bases: they become the tail clip, gap 1's D stays */
	in[0] = CIG_MK(64, OP_M);    /* contig [90, 154): ... 20 block2, 4 novel */
	rc = mm_lift_cigar(&R, 257, 90, in, 1, 0, 64, 0, &rso, &qso, &qeo, out, 32, &n_out);
	expect_int("ends in gap 2 rc", rc, 0);
	expect_cigar("ends in gap 2", out, n_out, "11M3I2D26M5D20M");
	expect_int("ends in gap 2 qe", qeo, 60);
	/* wholly inside gap 2's novel bases -> no representation */
	in[0] = CIG_MK(5, OP_M);
	rc = mm_lift_cigar(&R, 257, 151, in, 1, 0, 5, 0, &rso, &qso, &qeo, out, 32, &n_out);
	expect_int("wholly in gap 2 -> -2", rc, -2);
}

/* Randomised oracle.  The contig->base map is derived INDEPENDENTLY of
 * mm_lift: a sequential walk over (block length, ins, del) pieces fills an
 * explicit array map[c] = base coordinate or -1 (novel).  For random reads
 * (M runs with the read's own I/D) the lifted CIGAR is walked and every
 * query base checked: an M base must land on map[c] of the contig base it was
 * aligned to, a base aligned to a novel contig base must come out as I or as
 * clipping, and every D in the output must cover exactly base coordinates
 * that no anchored query base uses (a panel deletion or the read's own). */
static uint32_t rng_s = 12345;
static uint32_t rng(void) { rng_s = rng_s * 1103515245u + 12345u; return (rng_s >> 8) & 0xffffff; }

static void t12_oracle(void)
{
	int it, n_case = 0, n_bad = 0, n_multi_lifted = 0;
	puts("[12] randomised oracle: explicit per-base map vs mm_lift_cigar");
	for (it = 0; it < 20000; ++it) {
		mm_lift_gap_t gx[MM_LIFT_MAX_GAPS];
		mm_lift_rec_t rec;
		int32_t ng = 1 + rng() % MM_LIFT_MAX_GAPS, g, c, clen = 0, blen;
		int64_t map[2048], t, off_a = 100000 + rng() % 1000;
		int32_t jct[MM_LIFT_MAX_GAPS], ins[MM_LIFT_MAX_GAPS], del[MM_LIFT_MAX_GAPS];
		uint32_t in[64], out[128];
		int32_t n_in = 0, n_out, qso, qeo, qlen = 0, rs, span, k, q, rev = rng() & 1;
		int64_t rso;
		int rc;
		/* build the contig piece by piece, filling the explicit map */
		t = off_a;
		blen = 1 + rng() % 120;
		for (c = 0; c < blen; ++c) map[clen++] = t++;
		for (g = 0; g < ng; ++g) {
			int kind = rng() % 3;
			ins[g] = kind != 1 ? 1 + rng() % 30 : 0;
			del[g] = kind != 0 ? 1 + rng() % 30 : 0;
			jct[g] = clen;
			for (c = 0; c < ins[g]; ++c) map[clen++] = -1;
			t += del[g];
			blen = 1 + rng() % (g == ng - 1 ? 120 : 40);
			for (c = 0; c < blen; ++c) map[clen++] = t++;
		}
		for (g = 1; g < ng; ++g) { gx[g-1].jct = jct[g]; gx[g-1].ins_len = ins[g]; gx[g-1].del_len = del[g]; }
		memset(&rec, 0, sizeof rec);
		rec.base_rid = 0; rec.jct = jct[0]; rec.ins_len = ins[0]; rec.del_len = del[0]; rec.off_a = off_a;
		rec.n_xgap = ng - 1; rec.xgap = ng > 1 ? gx : 0;
		/* a random read: start, then M runs with occasional own I (1-3) / D (1-3) */
		rs = rng() % clen;
		span = 0;
		while (n_in < 60 && rs + span < clen) {
			int32_t m = 1 + rng() % 60;
			if (rs + span + m > clen) m = clen - rs - span;
			in[n_in++] = CIG_MK(m, OP_M); span += m; qlen += m;
			if (rng() % 4 == 0 || rs + span >= clen) break;
			if (rng() & 1) { int32_t l = 1 + rng() % 3; in[n_in++] = CIG_MK(l, OP_I); qlen += l; }
			else { int32_t l = 1 + rng() % 3; if (rs + span + l >= clen) break; in[n_in++] = CIG_MK(l, OP_D); span += l; }
		}
		while (n_in > 0 && CIG_OP(in[n_in-1]) != OP_M) { if (CIG_OP(in[n_in-1]) == OP_I) qlen -= CIG_LEN(in[n_in-1]); else span -= CIG_LEN(in[n_in-1]); --n_in; }
		if (n_in == 0) continue;
		++n_case;
		rc = mm_lift_cigar(&rec, clen, rs, in, n_in, 0, qlen, rev, &rso, &qso, &qeo, out, 128, &n_out);
		{
			/* expected: the forward-order list of (query index -> base coord or -1 novel) */
			int64_t qmap[2048]; int any = 0, ok = 1;
			int32_t lead, tail, first_anch = -1, last_anch = -1, oq, cc = rs;
			for (k = 0, q = 0; k < n_in; ++k) {
				int op = CIG_OP(in[k]); int32_t l = CIG_LEN(in[k]), j;
				if (op == OP_M) { for (j = 0; j < l; ++j) { qmap[q] = map[cc + j]; if (qmap[q] >= 0) { any = 1; if (first_anch < 0) first_anch = q; last_anch = q; } ++q; } cc += l; }
				else if (op == OP_I) { for (j = 0; j < l; ++j) qmap[q++] = -2; }   /* read's own insertion */
				else cc += l;
			}
			if (!any) { if (rc != -2) { ok = 0; } goto judge; }
			if (rc != 0) { ok = 0; goto judge; }
			/* clipping in CIGAR order = the query bases before the first / after the last anchored base */
			lead = first_anch; tail = qlen - 1 - last_anch;
			if (!rev) { if (qso != lead || qeo != qlen - tail) ok = 0; }
			else      { if (qso != tail || qeo != qlen - lead) ok = 0; }
			/* walk the output: M bases must match qmap exactly; I bases must be unanchored */
			{
				int64_t tt = rso, prev_t = -1; int32_t kk;
				oq = lead;
				for (kk = 0; kk < n_out && ok; ++kk) {
					int op = CIG_OP(out[kk]); int32_t l = CIG_LEN(out[kk]), j;
					if (l <= 0) ok = 0;
					if (kk > 0 && CIG_OP(out[kk-1]) == (uint32_t)op) ok = 0;            /* adjacent ops differ */
					if ((kk == 0 || kk == n_out - 1) && (op == OP_D || op == OP_I)) ok = 0;   /* no I/D at either end */
					if (op == OP_M) { for (j = 0; j < l; ++j) { if (qmap[oq] != tt) ok = 0; if (prev_t >= 0 && tt <= prev_t) ok = 0; prev_t = tt; ++oq; ++tt; } }
					else if (op == OP_I) { for (j = 0; j < l; ++j) { if (qmap[oq] >= 0) ok = 0; ++oq; } }
					else if (op == OP_D) tt += l;
					else ok = 0;
				}
				if (oq != qlen - tail) ok = 0;                                              /* query length preserved */
				if (tt != qmap[last_anch] + 1) ok = 0;                                      /* ends on the last anchored base */
			}
			if (ok && ng > 1) ++n_multi_lifted;
judge:
			if (!ok) {
				char b1[1024], b2[1024];
				if (++n_bad <= 5) {
					cig_str(in, n_in, b1); cig_str(out, rc == 0 ? n_out : 0, b2);
					printf("  FAIL oracle case %d: ng %d rs %d rev %d in %s -> rc %d out %s rs %lld q[%d,%d)\n", it, ng, rs, rev, b1, rc, b2, (long long)rso, qso, qeo);
				}
			}
		}
	}
	++checks;
	if (n_bad) ++fails;
	printf("  %d random cases, %d disagree with the explicit map; %d lifted across >1 gap\n", n_case, n_bad, n_multi_lifted);
	++checks;
	if (n_multi_lifted < 1000) { printf("  FAIL oracle exercised too few multi-gap cases (%d)\n", n_multi_lifted); ++fails; }
}

/* same locus: the single contig of A and the haplotype contig A+B, for a read
 * that carries only A and ends before B, project to the SAME alignment, and
 * for a read reaching past B (which on A+B shows B's gap undone by the read's
 * own ops) they still share block 0's diagonal -> one locus, not two. */
static void t13_same_locus_hap(void)
{
	uint32_t in[4], o1[32], o2[32];
	int32_t n1, n2, qs1, qe1, qs2, qe2;
	int64_t rs1, rs2;
	mm_reg1_t a, b; mm_extra_t *pa, *pb;
	puts("[13] same locus: single contig A vs haplotype contig A+B");
	in[0] = CIG_MK(89, OP_M);    /* read carries A, ends before B */
	mm_lift_cigar(&H2_A, 201, 51, in, 1, 0, 89, 0, &rs1, &qs1, &qe1, o1, 32, &n1);
	mm_lift_cigar(&H2,   245, 51, in, 1, 0, 89, 0, &rs2, &qs2, &qe2, o2, 32, &n2);
	expect_cigar("A-contig projection", o1, n1, "50M4I35M");
	expect_cigar("A+B-contig projection identical", o2, n2, "50M4I35M");
	expect_i64("same rs", rs1, rs2);
	a = mk_reg(0, 0, (int)rs1, qs1, qe1, o1, n1, &pa);
	b = mk_reg(0, 0, (int)rs2, qs2, qe2, o2, n2, &pb);
	expect_int("A vs A+B, read ends before B: same locus", mm_lift_same_locus(&a, &b), 1);
	free(pa); free(pb);
	/* read carries only A and runs past B: on A+B the aligner has to put B's
	 * 6 deleted bases back as the read's own insertion (94M 6I 50M) */
	in[0] = CIG_MK(150, OP_M);
	mm_lift_cigar(&H2_A, 205, 51, in, 1, 0, 150, 0, &rs1, &qs1, &qe1, o1, 32, &n1);
	in[0] = CIG_MK(94, OP_M); in[1] = CIG_MK(6, OP_I); in[2] = CIG_MK(50, OP_M);
	mm_lift_cigar(&H2, 245, 51, in, 3, 0, 150, 0, &rs2, &qs2, &qe2, o2, 32, &n2);
	expect_cigar("A-contig projection, past B", o1, n1, "50M4I96M");
	expect_cigar("A+B projection keeps both, the read undoes B", o2, n2, "50M4I40M6I6D50M");
	a = mk_reg(0, 0, (int)rs1, qs1, qe1, o1, n1, &pa);
	b = mk_reg(0, 0, (int)rs2, qs2, qe2, o2, n2, &pb);
	expect_int("A vs A+B, read past B: still one locus", mm_lift_same_locus(&a, &b), 1);
	free(pa); free(pb);
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
	t9_parse_hap();
	t10_two_gaps();
	t11_three_gaps();
	t12_oracle();
	t13_same_locus_hap();
	printf("\n%d checks, %d failures\n", checks, fails);
	return fails ? 1 : 0;
}
