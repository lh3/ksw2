#include <stdio.h> // for debugging only
#include "ksw2.h"

typedef struct { int32_t h, e; } eh_t;

#define NEG_INF -0x40000000

int ksw_gg(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t gapo, int8_t gape, int w, int *n_cigar_, uint32_t **cigar_)
{
	eh_t *eh;
	int8_t *qp; // query profile
	int32_t i, j, k, max_j = 0, gapoe = gapo + gape, score, n_col, *off = 0, last_en = -1;
	uint8_t *z = 0; // backtrack matrix; in each cell: f<<4|e<<2|h; in principle, we can halve the memory, but backtrack will be more complex

	// allocate memory
	n_col = qlen < 2*w+1? qlen : 2*w+1; // maximum #columns of the backtrack matrix
	qp = kmalloc(km, qlen * m);
	eh = kcalloc(km, qlen + 1, 8);
	if (n_cigar_ && cigar_) {
		*n_cigar_ = 0;
		z = kmalloc(km, (size_t)n_col * tlen);
		off = kcalloc(km, tlen, 4);
	}

	// generate the query profile
	for (k = i = 0; k < m; ++k) {
		const int8_t *p = &mat[k * m];
		for (j = 0; j < qlen; ++j) qp[i++] = p[query[j]];
	}

	// fill the first row
	eh[0].h = 0, eh[0].e = -gapoe - gapo;
	for (j = 1; j <= qlen && j <= w; ++j)
		eh[j].h = -(gapo + gape * j), eh[j].e = -(gapoe + gapo + gape * j);
	for (; j <= qlen; ++j) eh[j].h = eh[j].e = NEG_INF; // everything is -inf outside the band

	// DP loop
	for (i = 0; i < tlen; ++i) { // target sequence is in the outer loop
		int32_t f = NEG_INF, h1, st, en, max = NEG_INF;
		int8_t *q = &qp[target[i] * qlen];
		#if 0
		st = max_j > w? max_j - w : 0;
		en = max_j + w + 1 < qlen? max_j + w + 1 : qlen;
		#else
		st = i > w? i - w : 0;
		en = i + w + 1 < qlen? i + w + 1 : qlen;
		#endif
		h1 = st > 0? NEG_INF : -(gapo + gape * i);
		f  = st > 0? NEG_INF : -(gapoe + gapo + gape * i);
		off[i] = st;
		if (n_cigar_ && cigar_) {
			uint8_t *zi = &z[(long)i * n_col];
			for (j = st; j < en; ++j) {
				// At the beginning of the loop: eh[j] = { H(i-1,j-1), E(i,j) }, f = F(i,j) and h1 = H(i,j-1)
				// Cells are computed in the following order:
				//   H(i,j)   = max{H(i-1,j-1) + S(i,j), E(i,j), F(i,j)}
				//   E(i+1,j) = max{H(i,j)-gapo, E(i,j)} - gape
				//   F(i,j+1) = max{H(i,j)-gapo, F(i,j)} - gape
				eh_t *p = &eh[j];
				int32_t h = p->h, e = p->e;
				uint8_t d; // direction
				p->h = h1;
				h += q[j];
				d = h >= e? 0 : 1;
				h = h >= e? h : e;
				d = h >= f? d : 2;
				h = h >= f? h : f;
				h1 = h;
				max_j = max > h? max_j : j;
				max   = max > h? max   : h;
				h -= gapoe;
				e -= gape;
				d |= e > h? 1<<2 : 0;
				e  = e > h? e    : h;
				p->e = e;
				f -= gape;
				d |= f > h? 2<<4 : 0; // if we want to halve the memory, use one bit only, instead of two
				f  = f > h? f    : h;
				zi[j - st] = d; // z[i,j] keeps h for the current cell and e/f for the next cell
			}
		} else {
			for (j = st; j < en; ++j) {
				eh_t *p = &eh[j];
				int32_t h = p->h, e = p->e;
				p->h = h1;
				h += q[j];
				h = h >= e? h : e;
				h = h >= f? h : f;
				h1 = h;
				max_j = max > h? max_j : j;
				max   = max > h? max   : h;
				h -= gapoe;
				e -= gape;
				e  = e > h? e : h;
				p->e = e;
				f -= gape;
				f  = f > h? f : h;
			}
		}
		eh[en].h = h1, eh[en].e = NEG_INF, last_en = en;
	}

	// backtrack
	score = eh[qlen].h;
	if (n_cigar_ && cigar_) {
		int n_cigar = 0, m_cigar = 0, which = 0;
		uint32_t *cigar = 0, tmp;
		i = tlen - 1, k = last_en - 1; // (i,k) points to the last cell; FIXME: with a moving band, we need to take care of last deletion/insertion!!!
		while (i >= 0 && k >= 0) {
			tmp = z[i * n_col + k - off[i]];
			which = tmp >> (which << 1) & 3;
			if (which == 0 && tmp>>6) break;
			if (which == 0) which = tmp & 3;
			if (which == 0)      cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, 0, 1), --i, --k; // match
			else if (which == 1) cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, 2, 1), --i;      // deletion
			else                 cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, 1, 1), --k;      // insertion
		}
		if (i >= 0) cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, 2, i + 1); // first deletion
		if (k >= 0) cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, 1, k + 1); // first insertion
		for (i = 0; i < n_cigar>>1; ++i) // reverse CIGAR
			tmp = cigar[i], cigar[i] = cigar[n_cigar-1-i], cigar[n_cigar-1-i] = tmp;
		*n_cigar_ = n_cigar, *cigar_ = cigar;
	}

	kfree(km, qp); kfree(km, eh);
	if (n_cigar_ && cigar_) {
		kfree(km, z);
		kfree(km, off);
	}
	return score;
}
