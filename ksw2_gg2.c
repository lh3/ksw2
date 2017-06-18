#include <stdio.h> // for debugging only
#include "ksw2.h"

int ksw_gg2(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t q, int8_t e, int w, int *n_cigar_, uint32_t **cigar_)
{
	int qe = q + e, qe2 = qe + qe, r, t, n_col, *off, score;
	int8_t *u, *v, *x, *y, *s;
	uint8_t *p, *qr;

	u = (int8_t*)kcalloc(km, tlen + 1, 1);
	v = (int8_t*)kcalloc(km, tlen + 1, 1);
	x = (int8_t*)kcalloc(km, tlen + 1, 1);
	y = (int8_t*)kcalloc(km, tlen + 1, 1);
	s = (int8_t*)kmalloc(km, tlen);
	qr = (uint8_t*)kmalloc(km, qlen);
	n_col = w + 1 < tlen? w + 1 : tlen;
	p = (uint8_t*)kcalloc(km, (qlen + tlen) * n_col, 1);
	off = (int*)kmalloc(km, (qlen + tlen) * sizeof(int));

	for (t = 0; t < qlen; ++t)
		qr[t] = query[qlen - 1 - t];

	for (r = 0; r < qlen + tlen - 1; ++r) {
		int st = 0, en = tlen - 1;
		int8_t x1, v1;
		uint8_t *pr = p + r * n_col;
		// find the boundaries
		if (st < r - qlen + 1) st = r - qlen + 1;
		if (en > r) en = r;
		if (st < (r-w+1)>>1) st = (r-w+1)>>1; // take the ceil
		if (en > (r+w)>>1) en = (r+w)>>1; // take the floor
		off[r] = st;
		// set boundary conditions
		if (st != 0) {
			if (r > st + st + w - 1) x1 = v1 = 0;
			else x1 = x[st-1], v1 = v[st-1]; // (r-1, st-1) in the band
		} else x1 = 0, v1 = r? q : 0;
		if (en != r) {
			if (r < en + en - w - 1) y[en] = u[en] = 0; // (r-1,en) out of the band; TODO: is this line necessary?
		} else y[r] = 0, u[r] = r? q : 0;
		// loop fission: set scores first
		for (t = st; t <= en; ++t)
			s[t] = mat[target[t] * m + qr[t + qlen - 1 - r]];
		// core loop
		for (t = st; t <= en; ++t) {
			/* At the beginning of the loop, v1=v(r-1,t-1), x1=x(r-1,t-1), u[t]=u(r-1,t), v[t]=v(r-1,t), x[t]=x(r-1,t), y[t]=y(r-1,t)
			   a      = x(r-1,t-1) + v(r-1,t-1)
			   b      = y(r-1,t)   + u(r-1,t)
			   z      = max{ S(t,r-t) + 2q + 2r, a, b }
			   u(r,t) = z - v(r-1,t-1)
			   v(r,t) = z - u(r-1,t)
			   x(r,t) = max{ 0, a - z + q }
			   y(r,t) = max{ 0, b - z + q }
			*/
			uint8_t d;
			int8_t u1;
			int8_t z = s[t] + qe2;
			int8_t a = x1   + v1;
			int8_t b = y[t] + u[t];
			d = a > z? 1 : 0; // d = z >= a? 0 : 1
			z = a > z? a : z;
			d = b > z? 2 : d;
			z = b > z? b : z;
			u1 = u[t];              // u1   = u(r-1,t) (temporary variable)
			u[t] = z - v1;          // u[t] = u(r,t)
			v1 = v[t];              // v1   = v(r-1,t) (set for the next iteration)
			v[t] = z - u1;          // v[t] = v(r,t)
			z -= q;
			a -= z;
			b -= z;
			x1 = x[t];              // x1   = x(r-1,t) (set for the next iteration)
			d   |= a > 0? 1<<2 : 0;
			x[t] = a > 0? a    : 0; // x[t] = x(r,t)
			d   |= b > 0? 2<<4 : 0;
			y[t] = b > 0? b    : 0; // y[t] = y(r,t)
			pr[t - st] = d;
		}
	}
	kfree(km, u); kfree(km, v); kfree(km, x); kfree(km, y); kfree(km, s); kfree(km, qr);
	{ // backtrack
		int n_cigar = 0, m_cigar = 0, which = 0, i, j, k, l;
		uint32_t *cigar = 0, tmp;
		i = tlen - 1, j = qlen - 1;
		while (i >= 0 && j >= 0) {
			r = i + j;
			tmp = p[r * n_col + i - off[r]];
			which = tmp >> (which << 1) & 3;
			if (which == 0 && tmp>>6) break;
			if (which == 0) which = tmp & 3;
			if (which == 0)      cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, 0, 1), --i, --j; // match
			else if (which == 1) cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, 2, 1), --i;      // deletion
			else                 cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, 1, 1), --j;      // insertion
		}
		if (i >= 0) cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, 2, i + 1); // first deletion
		if (j >= 0) cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, 1, j + 1); // first insertion
		for (i = 0; i < n_cigar>>1; ++i) // reverse CIGAR
			tmp = cigar[i], cigar[i] = cigar[n_cigar-1-i], cigar[n_cigar-1-i] = tmp;
		*n_cigar_ = n_cigar, *cigar_ = cigar;

		// compute score
		for (k = 0, score = 0, i = j = 0; k < n_cigar; ++k) {
			int op = cigar[k] & 0xf, len = cigar[k] >> 4;
			if (op == 0) {
				for (l = 0; l < len; ++l)
					score += mat[target[i + l] * m + query[j + l]];
				i += len, j += len;
			} else if (op == 1) {
				score -= q + len * e;
				j += len;
			} else if (op == 2) {
				score -= q + len * e;
				i += len;
			}
		}
	}
	kfree(km, p); kfree(km, off);
	return score;
}
