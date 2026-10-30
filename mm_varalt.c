/* --var-alt: known-ALT-aware re-scoring of candidates (see mm_varalt.h). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "mmpriv.h"
#include "mm_varalt.h"

static const int8_t va_nt4[256] = {
	['A'] = 1, ['C'] = 2, ['G'] = 3, ['T'] = 4, ['a'] = 1, ['c'] = 2, ['g'] = 3, ['t'] = 4
}; /* code + 1; 0 = not ACGT */

void mm_varalt_destroy(mm_varalt_t *va)
{
	if (va == NULL) return;
	free(va->s);
	free(va->pos_pool);
	free(va->alt_pool);
	free(va);
}

mm_varalt_t *mm_varalt_load(const mm_idx_t *mi, const char *fn)
{
	gzFile fp;
	mm_varalt_t *va;
	char *line;
	int line_cap = 1 << 16, last_rid = -1;
	char last_name[1024];
	uint64_t *key = NULL;
	int64_t n_key = 0, m_key = 0, i, k;
	uint32_t *cnt;

	if (mi->h == 0) mm_idx_index_name((mm_idx_t*)mi);
	fp = strcmp(fn, "-") == 0 ? gzdopen(0, "r") : gzopen(fn, "r");
	if (fp == 0) return NULL;
	va = (mm_varalt_t*)calloc(1, sizeof(mm_varalt_t));
	va->n_seq = mi->n_seq;
	va->s = (mm_va_seq_t*)calloc(mi->n_seq, sizeof(mm_va_seq_t));
	line = (char*)malloc(line_cap);
	last_name[0] = 0;
	while (gzgets(fp, line, line_cap) != NULL) {
		char *p = line, *name, *end;
		long pos;
		int rid, mask = 0, l;
		l = strlen(line);
		if (l == line_cap - 1 && line[l-1] != '\n') { /* overlong line: skip */
			++va->n_bad;
			while (gzgets(fp, line, line_cap) != NULL) { l = strlen(line); if (l < line_cap - 1 || line[l-1] == '\n') break; }
			continue;
		}
		while (l > 0 && (line[l-1] == '\n' || line[l-1] == '\r')) line[--l] = 0;
		if (l == 0 || line[0] == '#') continue;
		++va->n_lines;
		name = p;
		while (*p && *p != '\t') ++p;
		if (*p != '\t') { ++va->n_bad; continue; }
		*p++ = 0;
		pos = strtol(p, &end, 10);
		if (end == p || *end != '\t' || pos < 1) { ++va->n_bad; continue; }
		p = end + 1;
		for (; *p && *p != '\t'; ++p) {
			int c = va_nt4[(uint8_t)*p];
			if (c == 0) { mask = -1; break; }
			mask |= 1 << (c - 1);
		}
		if (mask <= 0) { ++va->n_bad; continue; }
		if (last_rid >= 0 && strcmp(name, last_name) == 0) rid = last_rid;
		else {
			rid = mm_idx_name2id(mi, name);
			if (strlen(name) < sizeof(last_name)) { strcpy(last_name, name); last_rid = rid; }
			else last_name[0] = 0, last_rid = -1;
		}
		if (rid < 0) { ++va->n_unknown_seq; continue; }
		if ((uint64_t)pos > mi->seq[rid].len) { ++va->n_out_of_range; continue; }
		{	/* drop an ALT equal to the base */
			uint8_t b;
			mm_idx_getseq(mi, rid, pos - 1, pos, &b);
			if (b < 4 && (mask >> b & 1)) { ++va->n_eq_ref; mask &= ~(1 << b); }
			if (mask == 0) continue;
		}
		if (n_key == m_key) {
			m_key = m_key ? m_key + (m_key >> 1) : 1 << 20;
			key = (uint64_t*)realloc(key, m_key * sizeof(uint64_t));
		}
		key[n_key++] = (uint64_t)rid << 36 | (uint64_t)(pos - 1) << 4 | (uint64_t)mask;
	}
	gzclose(fp);
	free(line);
	radix_sort_64(key, key + n_key);
	/* merge equal (rid, pos), count per sequence */
	for (i = 1, k = 0; i < n_key; ++i) {
		if (key[i] >> 4 == key[k] >> 4) { key[k] |= key[i] & 0xf; ++va->n_dup; }
		else key[++k] = key[i];
	}
	n_key = n_key ? k + 1 : 0;
	va->n_pos = n_key;
	va->pos_pool = (uint32_t*)malloc((n_key ? n_key : 1) * sizeof(uint32_t));
	va->alt_pool = (uint8_t*)malloc(n_key ? n_key : 1);
	cnt = (uint32_t*)calloc(mi->n_seq, sizeof(uint32_t));
	for (i = 0; i < n_key; ++i) {
		int m = key[i] & 0xf;
		++cnt[key[i] >> 36];
		va->pos_pool[i] = (uint32_t)(key[i] >> 4);   /* truncation keeps the low 32 bits = pos */
		va->alt_pool[i] = m;
		va->n_alleles += (m & 1) + (m >> 1 & 1) + (m >> 2 & 1) + (m >> 3 & 1);
	}
	for (i = 0, k = 0; i < (int64_t)mi->n_seq; ++i) {
		va->s[i].n = cnt[i];
		va->s[i].pos = va->pos_pool + k;
		va->s[i].alt = va->alt_pool + k;
		k += cnt[i];
	}
	free(cnt);
	free(key);
	return va;
}

static inline uint32_t va_lower_bound(const mm_va_seq_t *s, uint32_t x)
{
	uint32_t lo = 0, hi = s->n;
	while (lo < hi) {
		uint32_t mid = lo + ((hi - lo) >> 1);
		if (s->pos[mid] < x) lo = mid + 1;
		else hi = mid;
	}
	return lo;
}

int mm_varalt_mask(const mm_varalt_t *va, int32_t rid, uint32_t pos)
{
	const mm_va_seq_t *s;
	uint32_t j;
	if (va == NULL || rid < 0 || (uint32_t)rid >= va->n_seq) return 0;
	s = &va->s[rid];
	j = va_lower_bound(s, pos);
	return j < s->n && s->pos[j] == pos ? s->alt[j] : 0;
}

int mm_varalt_rescore(const mm_idx_t *mi, const mm_mapopt_t *opt, int qlen, const char *qstr, mm_reg1_t *r)
{
	extern unsigned char seq_nt4_table[256];
	const mm_varalt_t *va = mi->varalt;
	const mm_va_seq_t *s;
	mm_extra_t *p = r->p;
	uint8_t *qbuf, *qseq, *tseq;
	uint32_t j, k, l;
	int32_t i, l_q, l_t, qoff = 0, toff = 0, n = 0, amb;
	double s0 = 0.0, s1 = 0.0, max0 = 0.0, max1 = 0.0;
	int log_gap = !(opt->flag & (MM_F_SR|MM_F_SR_RNA));   /* as mm_align1() calls mm_update_extra() */

	(void)qlen;
	r->va = 0;
	if (va == NULL || p == NULL || p->n_cigar == 0 || r->inv) return 0;
	if (r->rid < 0 || (uint32_t)r->rid >= va->n_seq) return 0;
	s = &va->s[r->rid];
	if (s->n == 0) return 0;
	j = va_lower_bound(s, (uint32_t)r->rs);
	if (j == s->n || s->pos[j] >= (uint32_t)r->re) return 0;   /* no known SNV under the alignment */
	l_q = r->qe - r->qs, l_t = r->re - r->rs;
	if (l_q <= 0 || l_t <= 0) return 0;
	qbuf = (uint8_t*)malloc((size_t)l_q + (size_t)l_t);
	if (qbuf == NULL) return 0;
	qseq = qbuf, tseq = qbuf + l_q;
	/* query in alignment orientation; r->qs/qe are forward-read coordinates */
	if (!r->rev) {
		for (i = 0; i < l_q; ++i) qseq[i] = seq_nt4_table[(uint8_t)qstr[r->qs + i]];
	} else {
		for (i = 0; i < l_q; ++i) {
			uint8_t c = seq_nt4_table[(uint8_t)qstr[r->qe - 1 - i]];
			qseq[i] = c < 4 ? 3 - c : 4;
		}
	}
	mm_idx_getseq(mi, r->rid, r->rs, r->re, tseq);
	amb = opt->sc_ambi > 0 ? -opt->sc_ambi : opt->sc_ambi;
	for (k = 0; k < p->n_cigar; ++k) {
		uint32_t op = p->cigar[k] & 0xf, len = p->cigar[k] >> 4;
		if (op == MM_CIGAR_MATCH || op == MM_CIGAR_EQ_MATCH || op == MM_CIGAR_X_MISMATCH) {
			for (l = 0; l < len; ++l) {
				int cq = qseq[qoff + l], ct = tseq[toff + l], sc0, sc1;
				if (cq > 3 || ct > 3) sc0 = sc1 = amb;
				else if (cq == ct) sc0 = sc1 = opt->a;
				else {
					uint32_t tp = (uint32_t)(r->rs + toff + l);
					sc0 = sc1 = -opt->b;
					while (j < s->n && s->pos[j] < tp) ++j;
					if (j < s->n && s->pos[j] == tp && (s->alt[j] >> cq & 1)) sc1 = opt->a, ++n;
				}
				s0 += sc0, s1 += sc1;
				if (s0 < 0) s0 = 0; else if (s0 > max0) max0 = s0;
				if (s1 < 0) s1 = 0; else if (s1 > max1) max1 = s1;
			}
			qoff += len, toff += len;
		} else if (op == MM_CIGAR_INS || op == MM_CIGAR_DEL) {
			double g = log_gap ? opt->q + (double)opt->e * mg_log2(1.0 + len) : opt->q + opt->e;
			s0 -= g, s1 -= g;
			if (s0 < 0) s0 = 0;
			if (s1 < 0) s1 = 0;
			if (op == MM_CIGAR_INS) qoff += len; else toff += len;
		} else if (op == MM_CIGAR_N_SKIP) toff += len;
	}
	free(qbuf);
	if (n == 0) return 0;
	p->dp_score += n * (opt->a + opt->b);
	p->dp_max += (int32_t)(max1 + .499) - (int32_t)(max0 + .499);
	r->mlen += n;
	r->va = n;
	return n;
}
