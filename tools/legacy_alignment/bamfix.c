/* Split a raw BAM into records whose CIGAR query length matches l_seq (OUT_OK) and
 * broken ones (OUT_BAD), which are rewritten as unmapped: CIGAR dropped, tid/pos = -1,
 * flag |= 0x4, proper-pair bit cleared. Little-endian hosts only; the caller
 * validates BAM framing beforehand. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "htslib/bgzf.h"
static void die(const char *m){ fprintf(stderr, "bamfix: %s\n", m); exit(1); }
static void copy_header(BGZF *in, BGZF *o1, BGZF *o2){
  char magic[4]; int32_t l_text, n_ref, i;
  if (bgzf_read(in, magic, 4) != 4 || memcmp(magic, "BAM\1", 4)) die("not a BAM file");
  bgzf_write(o1, magic, 4); bgzf_write(o2, magic, 4);
  bgzf_read(in, &l_text, 4); bgzf_write(o1, &l_text, 4); bgzf_write(o2, &l_text, 4);
  char *t = malloc(l_text); bgzf_read(in, t, l_text); bgzf_write(o1, t, l_text); bgzf_write(o2, t, l_text); free(t);
  bgzf_read(in, &n_ref, 4); bgzf_write(o1, &n_ref, 4); bgzf_write(o2, &n_ref, 4);
  for (i = 0; i < n_ref; i++) { int32_t l; bgzf_read(in, &l, 4); bgzf_write(o1, &l, 4); bgzf_write(o2, &l, 4);
    char *n = malloc(l + 4); bgzf_read(in, n, l + 4); bgzf_write(o1, n, l + 4); bgzf_write(o2, n, l + 4); free(n); }
}
int main(int argc, char **argv){
  if (argc != 4) die("bamfix IN.bam OUT_OK.bam OUT_BAD.bam");
  BGZF *in = bgzf_open(argv[1], "r"), *ok = bgzf_open(argv[2], "w"), *bad = bgzf_open(argv[3], "w");
  if (!in || !ok || !bad) die("cannot open files");
  bgzf_mt(ok, 4, 256); bgzf_mt(bad, 2, 256);
  copy_header(in, ok, bad);
  uint8_t *buf = malloc(1 << 20); size_t cap = 1 << 20; long n = 0, nbad = 0; int32_t bs;
  while (bgzf_read(in, &bs, 4) == 4) {
    if ((size_t)bs > cap) { cap = bs * 2; buf = realloc(buf, cap); }
    if (bgzf_read(in, buf, bs) != bs) die("truncated record");
    n++;
    uint8_t lrn = buf[8]; uint16_t ncig; memcpy(&ncig, buf + 12, 2); int32_t lseq; memcpy(&lseq, buf + 16, 4);
    uint32_t *cig = (uint32_t *)(buf + 32 + lrn); int64_t ql = 0; int i;
    for (i = 0; i < ncig; i++) { int op = cig[i] & 0xf; uint32_t len = cig[i] >> 4; if (op == 0 || op == 1 || op == 4 || op == 7 || op == 8) ql += len; }
    if (ncig == 0 || ql == lseq) { bgzf_write(ok, &bs, 4); bgzf_write(ok, buf, bs); continue; }
    nbad++;
    /* Unmapped copy: CIGAR removed, MAPQ and n_cigar zeroed, bin 4680 = reg2bin(-1, 0). */
    int32_t nbs = bs - 4 * ncig; uint8_t *o = malloc(nbs);
    int32_t m1 = -1; memcpy(o, &m1, 4); memcpy(o + 4, &m1, 4);          
    o[8] = lrn; o[9] = 0; uint16_t bin = 4680; memcpy(o + 10, &bin, 2); uint16_t z = 0; memcpy(o + 12, &z, 2);
    uint16_t fl; memcpy(&fl, buf + 14, 2); fl |= 4; fl &= ~0x2; memcpy(o + 14, &fl, 2);
    memcpy(o + 16, buf + 16, 16);                                          
    memcpy(o + 32, buf + 32, lrn);                                         
    memcpy(o + 32 + lrn, buf + 32 + lrn + 4 * ncig, bs - 32 - lrn - 4 * ncig);  
    bgzf_write(bad, &nbs, 4); bgzf_write(bad, o, nbs); free(o);
  }
  bgzf_close(in); bgzf_close(ok); bgzf_close(bad);
  printf("bamfix: %ld records, %ld with CIGAR/l_seq mismatch written as unmapped to %s\n", n, nbad, argv[3]);
  return 0;
}
