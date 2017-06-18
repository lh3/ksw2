#include <stdio.h> // for debugging only
#include "ksw2.h"

#ifdef __SSE2__
#include <emmintrin.h>

#ifdef __SSE4_1__
#include <smmintrin.h>
#endif

int ksw_gg2_sse(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t q, int8_t e, int w, int *n_cigar_, uint32_t **cigar_)
{
	int r, t, n_col, *off, score, tlen16;
	int8_t *u, *v, *x, *y, *s;
	uint8_t *p, *qr, *mem;
	__m128i q_, qe2_, zero_, flag1_, flag2_, flag4_, flag32_;

	zero_   = _mm_set1_epi8(0);
	q_      = _mm_set1_epi8(q);
	qe2_    = _mm_set1_epi8((q + e) * 2);
	flag1_  = _mm_set1_epi8(1<<0);
	flag2_  = _mm_set1_epi8(2<<0);
	flag4_  = _mm_set1_epi8(1<<2);
	flag32_ = _mm_set1_epi8(2<<4);

	w = (w + 1 + 15) / 16 * 16 - 1;
	tlen16 = (tlen + 15) / 16 * 16;
	n_col = w + 1 < tlen16? w + 1 : tlen16; // number of columns in the backtrack matrix
	n_col += 16, tlen16 += 16; // leave enough space at the end

	mem = (uint8_t*)kcalloc(km, tlen16 * 5 + 15, 1);
	u = (int8_t*)(((size_t)mem + 15) >> 4 << 4); // 16-byte aligned (though not necessary)
	v = u + tlen16, x = v + tlen16, y = x + tlen16, s = y + tlen16;
	qr = (uint8_t*)kcalloc(km, qlen, 1);
	p = (uint8_t*)kcalloc(km, (qlen + tlen) * n_col, 1);
	off = (int*)kmalloc(km, (qlen + tlen) * sizeof(int));

	for (t = 0; t < qlen; ++t)
		qr[t] = query[qlen - 1 - t];

	for (r = 0; r < qlen + tlen - 1; ++r) {
		int st = 0, en = tlen - 1;
		int8_t x1, v1;
		uint8_t *pr = p + r * n_col;
		__m128i x1_, v1_;
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
		x1_ = _mm_cvtsi32_si128(x1);
		v1_ = _mm_cvtsi32_si128(v1);
		for (t = st; t <= en; t += 16) {
			__m128i d, z, a, b, xt1, vt1, ut, tmp;

			z = _mm_add_epi8(_mm_loadu_si128((__m128i*)&s[t]), qe2_);

			xt1 = _mm_loadu_si128((__m128i*)&x[t]);          // xt1 <- x[r-1][t..t+15]
			tmp = _mm_srli_si128(xt1, 15);                   // tmp <- x[r-1][t+15]
			xt1 = _mm_or_si128(_mm_slli_si128(xt1, 1), x1_); // xt1 <- x[r-1][t-1..t+14]
			x1_ = tmp;
			vt1 = _mm_loadu_si128((__m128i*)&v[t]);          // vt1 <- v[r-1][t..t+15]
			tmp = _mm_srli_si128(vt1, 15);                   // tmp <- v[r-1][t+15]
			vt1 = _mm_or_si128(_mm_slli_si128(vt1, 1), v1_); // vt1 <- v[r-1][t-1..t+14]
			v1_ = tmp;
			a = _mm_add_epi8(xt1, vt1);                      // a <- x[r-1][t-1..t+14] + v[r-1][t-1..t+14]

			ut = _mm_loadu_si128((__m128i*)&u[t]);           // ut <- u[t..t+15]
			b = _mm_add_epi8(_mm_loadu_si128((__m128i*)&y[t]), ut); // b <- y[r-1][t..t+15] + u[r-1][t..t+15]

			d = _mm_and_si128(_mm_cmpgt_epi8(a, z), flag1_); // d = a > z? 1 : 0
#ifdef __SSE4_1__
			z = _mm_max_epi8(z, a);                          // z = z > a? z : a (signed)
			tmp = _mm_cmpgt_epi8(b, z);
			d = _mm_blendv_epi8(d, flag2_, tmp);             // d = b > z? 2 : d
#else // we need to emulate SSE4.1 intrinsics _mm_max_epi8() and _mm_blendv_epi8()
			z = _mm_and_si128(z, _mm_cmpgt_epi8(z, zero_));  // z = z > 0? z : 0;
			z = _mm_max_epu8(z, a);                          // z = max(z, a); this works because both are non-negative
			tmp = _mm_cmpgt_epi8(b, z);
			d = _mm_or_si128(_mm_andnot_si128(tmp, d), _mm_and_si128(tmp, flag2_)); // d = b > z? 2 : d; emulating blendv
#endif
			z = _mm_max_epu8(z, b);                          // z = max(z, b); this works because both are non-negative
			_mm_storeu_si128((__m128i*)&u[t], _mm_sub_epi8(z, vt1)); // u[r][t..t+15] <- z - v[r-1][t-1..t+14]
			_mm_storeu_si128((__m128i*)&v[t], _mm_sub_epi8(z, ut));  // v[r][t..t+15] <- z - u[r-1][t..t+15]

			z = _mm_sub_epi8(z, q_);
			a = _mm_sub_epi8(a, z);
			b = _mm_sub_epi8(b, z);
			tmp = _mm_cmpgt_epi8(a, zero_);
			d = _mm_or_si128(d, _mm_and_si128(flag4_,  tmp));
			_mm_storeu_si128((__m128i*)&x[t], _mm_and_si128(a, tmp));
			tmp = _mm_cmpgt_epi8(b, zero_);
			d = _mm_or_si128(d, _mm_and_si128(flag32_, tmp));
			_mm_storeu_si128((__m128i*)&y[t], _mm_and_si128(b, tmp));
			_mm_storeu_si128((__m128i*)&pr[t - st], d);
		}
		// for (t = st; t <= en; ++t) printf("(%d,%d)\t(%d,%d,%d,%d)\t%x\n", r, t, u[t], v[t], x[t], y[t], pr[t-st]); // for debugging
	}
	kfree(km, mem); kfree(km, qr);
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
#endif // __SSE2__
