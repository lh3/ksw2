#include "ksw2.h"

typedef struct { int32_t h, e; } eh_t;

int ksw_ext(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int gapo, int gape, int w,
			int end_bonus, int zdrop, int h0, int *_qle, int *_tle, int *_gtle, int *_gscore)
{
	eh_t *eh; // score array
	int8_t *qp; // query profile
	int i, j, k, gapoe = gapo + gape, st, en, max, max_i, max_j, max_j0, max_gap, max_ie, gscore;

	// allocate memory
	qp = kmalloc(km, (long)qlen * m);
	eh = kcalloc(km, qlen + 1, 8);

	// generate the query profile
	for (k = i = 0; k < m; ++k) {
		const int8_t *p = &mat[k * m];
		for (j = 0; j < qlen; ++j) qp[i++] = p[query[j]];
	}

	// fill the first row
	eh[0].h = h0, eh[0].e = h0 - gapoe - gapoe > 0? h0 - gapoe - gapoe : 0;
	for (j = 1; j <= qlen && j <= w; ++j) {
		eh[j].h = h0 - (gapoe + gape * j), eh[j].e = h0 - (gapoe + gapoe + gape * j);
		if (eh[j].e < 0) eh[j].e = 0;
		if (eh[j].h < 0) {
			eh[j].h = 0;
			break;
		}
	}

	// adjust $w if it is too large
	k = m * m;
	for (i = 0, max = 0; i < k; ++i) // get the max score
		max = max > mat[i]? max : mat[i];
	max_gap = (int)((double)((qlen < tlen? qlen : tlen) * max + end_bonus - gapo) / gape + 1.);
	max_gap = max_gap > 1? max_gap : 1;
	w = w < max_gap? w : max_gap;

	// DP loop
	max = h0, max_i = max_j = -1; max_ie = -1, gscore = -1; max_j0 = 0;
	st = 0, en = qlen;
	for (i = 0; i < tlen; ++i) {
		int f = 0, h1, m0 = 0;
		int8_t *q = &qp[target[i] * qlen];
		// apply the band and the constraint (if provided)
		if (st < max_j0 - w)     st = max_j0 - w;
		if (en > max_j0 + w + 1) en = max_j0 + w + 1;
		// compute the first column
		if (st == 0) {
			h1 = h0 - (gapoe + gape * i);
			if (h1 < 0) h1 = 0;
			f = h0 - (gapoe + gapoe + gape * i);
			if (f < 0) f = 0;
		} else h1 = 0;
		for (j = st; j < en; ++j) {
			// At the beginning of the loop: eh[j] = { H(i-1,j-1), E(i,j) }, f = F(i,j) and h1 = H(i,j-1)
			// Similar to SSE2-SW, cells are computed in the following order:
			//   H(i,j)   = max{H(i-1,j-1)+S(i,j), E(i,j), F(i,j)}
			//   E(i+1,j) = max{H(i,j)-gapo, E(i,j)} - gape
			//   F(i,j+1) = max{H(i,j)-gapo, F(i,j)} - gape
			eh_t *p = &eh[j];
			int h = p->h, e = p->e; // get H(i-1,j-1) and E(i-1,j)
			p->h = h1;          // set H(i,j-1) for the next row
			h += q[j];
			h = h > e? h : e;   // e and f are guaranteed to be non-negative, so h>=0 even if h<0
			h = h > f? h : f;
			h1 = h;             // save H(i,j) to h1 for the next column
			max_j0 = m0 > h? max_j0 : j; // record the position where max score is achieved
			m0     = m0 > h? m0     : h; // m0 is stored at eh[mj+1]
			h -= gapoe;
			h = h > 0? h : 0;
			e -= gape;
			e = e > h? e : h;   // computed E(i+1,j)
			p->e = e;           // save E(i+1,j) for the next row
			f -= gape;
			f = f > h? f : h;   // computed F(i,j+1)
		}
		eh[en].h = h1; eh[en].e = 0;
		if (j == qlen) {
			max_ie = gscore > h1? max_ie : i;
			gscore = gscore > h1? gscore : h1;
		}
		if (m0 == 0) break;
		if (m0 > max) {
			max = m0, max_i = i, max_j = max_j0;
		} else if (zdrop > 0) {
			int diff = (i - max_i) - (max_j0 - max_j);
			if (max - m0 - (diff < 0? -diff : diff) * gape > zdrop) break;
		}
		// update beg and end for the next round
		for (j = st; j < en && eh[j].h == 0 && eh[j].e == 0; ++j);
		st = j;
		for (j = en; j >= st && eh[j].h == 0 && eh[j].e == 0; --j);
		en = j + 2 < qlen? j + 2 : qlen;
		//beg = 0; end = qlen; // uncomment this line for debugging
	}
	kfree(km, eh); kfree(km, qp);
	if (_qle) *_qle = max_j + 1;
	if (_tle) *_tle = max_i + 1;
	if (_gtle) *_gtle = max_ie + 1;
	if (_gscore) *_gscore = gscore;
	return max;
}
