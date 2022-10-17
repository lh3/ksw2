// Created by Haoyang
//First convert this function to C
// and prepare for CUDA parallelization

#include <string.h>
#include <stdio.h>
#include <assert.h>
#include "ksw2.h"

#ifdef __SSE2__
#ifdef USE_SIMDE
#include <simde/x86/sse2.h>
#else
#include <emmintrin.h>
#endif

#ifdef KSW_SSE2_ONLY
#undef __SSE4_1__
#endif

#ifdef __SSE4_1__
#ifdef USE_SIMDE
#include <simde/x86/sse4.1.h>
#else
#include <smmintrin.h>
#endif
#endif


void ksw_extd2_c(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
				   int8_t q, int8_t e, int8_t q2, int8_t e2, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez)
{

	int r, t, qe = q + e, n_col, *off = 0, *off_end = 0, tlen_, qlen_, last_st, last_en, wl, wr, max_sc, min_sc, long_thres, long_diff;
	int with_cigar = !(flag&KSW_EZ_SCORE_ONLY), approx_max = !!(flag&KSW_EZ_APPROX_MAX);
	int32_t *H = 0, H0 = 0, last_H0_t = 0;
	uint8_t *qr, *sf, *mem, *mem2 = 0;
	int8_t q_, q2_, qe_, qe2_, zero_, sc_mch_, sc_mis_, m1_, sc_N_;
	int8_t *u, *v, *x, *y, *x2, *y2;
	int8_t *s;
	uint8_t *p = 0;

	ksw_reset_extz(ez);
	if (m <= 1 || qlen <= 0 || tlen <= 0) return;

	if (q2 + e2 < q + e) t = q, q = q2, q2 = t, t = e, e = e2, e2 = t; // make sure q+e no larger than q2+e2

	zero_   = 0;
	q_      = q;
	q2_     = q2;
	qe_     = q + e;
	qe2_    = q2 + e2;
	sc_mch_ = mat[0];
	sc_mis_ = mat[1];
	sc_N_   = mat[m*m-1] == 0? -e2 : mat[m*m-1];
	m1_     = m - 1; // wildcard

	if (w < 0) w = tlen > qlen? tlen : qlen;
	wl = wr = w;
	//tlen_ = (tlen + 15) / 16;    // tlen_: need how many vectors for reference chain
	n_col = qlen < tlen? qlen : tlen;  
	n_col = n_col < (w+1) ? n_col : w+1;
	for (t = 1, max_sc = mat[0], min_sc = mat[1]; t < m * m; ++t) {
		max_sc = max_sc > mat[t]? max_sc : mat[t];
		min_sc = min_sc < mat[t]? min_sc : mat[t];
	}   // process the mat, get max_sc min_sc
	if (-min_sc > 2 * (q + e)) return; // otherwise, we won't see any mismatches

	long_thres = e != e2? (q2 - q) / (e - e2) - 1 : 0;
	if (q2 + e2 + long_thres * e2 > q + e + long_thres * e)
		++long_thres;
	long_diff = long_thres * (e - e2) - (q2 - q) - e2;

	mem = (uint8_t*)kmalloc(km, tlen * 8 + qlen + 1); //?
	u = mem; // 16-byte aligned
	v = u + tlen, x = v + tlen, y = x + tlen, x2 = y + tlen, y2 = x2 + tlen;
	s = y2 + tlen, sf = s + tlen, qr = sf + tlen; //u,v,x,y,x2,y2,s,sf each length is tlen_*16 bytes, placed in mem
	memset(u,  -q  - e,  tlen);
	memset(v,  -q  - e,  tlen);
	memset(x,  -q  - e,  tlen);
	memset(y,  -q  - e,  tlen);
	memset(x2, -q2 - e2, tlen);
	memset(y2, -q2 - e2, tlen);
	if (!approx_max) {
		H = (int32_t*)kmalloc(km, tlen * 4);
		for (t = 0; t < tlen; ++t) H[t] = KSW_NEG_INF;  //Init H[i] = -inf
	}
	if (with_cigar) {
		mem2 = (uint8_t*)kmalloc(km, ((size_t)(qlen + tlen - 1) * n_col + 1)); // qlen + tlen - 1 = number of diagonal iterations
		p = (uint8_t*)mem2; //16-byte aligned, represent using p
		// In the backtrack matrix, value p[] has the following structure:
		//   bit 0-2: which type gets the max - 0 for H, 1 for E, 2 for F, 3 for \tilde{E} and 4 for \tilde{F}
		//   bit 3/0x08: 1 if a continuation on the E state (bit 5/0x20 for a continuation on \tilde{E})
		//   bit 4/0x10: 1 if a continuation on the F state (bit 6/0x40 for a continuation on \tilde{F})
		off = (int*)kmalloc(km, (qlen + tlen - 1) * sizeof(int) * 2); // off is a int array using r to index
		off_end = off + qlen + tlen - 1;
	}

	for (t = 0; t < qlen; ++t) qr[t] = query[qlen - 1 - t]; //qr is query reverse
	memcpy(sf, target, tlen); // sf = reference chain

	for (r = 0, last_st = last_en = -1; r < qlen + tlen - 1; ++r) {   //Outer for loop, last_st = "last start", "last end"
	//for (r = 0, last_st = last_en = -1; r < 0; ++r) {   // debug
		int st = 0, en = tlen - 1, st0, en0; // , st_, en_; // r is iteration index of organal; t_len is at x axis, q_len is the atg y_axis
		int8_t x1, x21, v1;
		uint8_t *qrr = qr + (qlen - 1 - r);  // This is like a offset pointer. Effect: qrr[t] & sf[t] is the needed index when we have  st < t < en.
		int8_t *u8 = (int8_t*)u, *v8 = (int8_t*)v;  //, *x8 = (int8_t*)x, *x28 = (int8_t*)x2;  //x8 is 8-byte address of x
		int8_t x1_, x21_, v1_;
		// find the boundaries
		if (st < r - qlen + 1) st = r - qlen + 1;  // Correct the value of st iif r is large, during each iteration (Can be explained using figure)
		if (en > r) en = r; // Correct the en if r is small
		if (st < (r-wr+1)>>1) st = (r-wr+1)>>1; // take the ceil, choose the band position
		if (en > (r+wl)>>1) en = (r+wl)>>1; // take the floor, choose the band position
		if (st > en) {
			ez->zdropped = 1;
			break;
		}
		st0 = st, en0 = en;   // Store the orignal st and en into st0 and en0

		if (st > 0) {
			if (st - 1 >= last_st && st - 1 <= last_en) {
				x1 = x[st - 1], x21 = x2[st - 1], v1 = v[st - 1]; // (r-1,s-1) calculated in the last round
			} else {  // else means cannot use at least one "pre-value" in the last round
				x1 = -q - e, x21 = -q2 - e2;
				v1 = -q - e;   // assign initial value
			}
		} else {  //st == 0
			x1 = -q - e, x21 = -q2 - e2;
			v1 = r == 0? -q - e : r < long_thres? -e : r == long_thres? long_diff : -e2;
		}
		if (en >= r) {  // en==r
			y[r] = -q - e, y2[r] = -q2 - e2;
			u[r] = r == 0? -q - e : r < long_thres? -e : r == long_thres? long_diff : -e2;
		}
		// loop fission: set scores first
		if (!(flag & KSW_EZ_GENERIC_SC)) {
			for (t = st; t <= en; t ++) {
				int8_t sq, st, score;
				sq = sf[t]; //sq = sf = reference chain
				st = qrr[t]; //st = qrr = partial reverse query chain
				if (sq==m1_||st==m1_)
				{
					score = sc_N_;
				}
				else{
					if (sq==st)
					{
						score = sc_mch_;
					}
					else
					{
						score = sc_mis_;
					}
				}
				*(s+t) = score;
			}

		} else {
			for (t = st; t <= en; ++t)
				((uint8_t*)s)[t] = mat[sf[t] * m + qrr[t]];  // Preprocess the score, matrix is like a lookup table.
		}
		// core loop
		x1_  = x1;
		x21_ = x21;
		v1_  = v1;
		assert(en - st + 1 <= n_col);

		// SCORE_ONLY
		// 0x01
		if (!with_cigar) { // score only
			int8_t u_new, v_new, x_new, y_new, x2_new, y2_new;
			int8_t u_next, v_next, x_next, y_next, x2_next, y2_next;
			for (t = st; t <= en; ++t) {
				int8_t z, a, b, a2, b2, xt1, x2t1, vt1, ut, tmp;

				z = s[t];
				xt1 = t == st ? x1_ : x[t-1];
				x1 = x[t];
				vt1 = t == st ? v1_ : v[t-1];
				v1 = v[t];
				a = xt1 + vt1;
				ut = u[t];
				b = y[t] + u[t];
				x2t1 = t == st ? x21_ : x2[t-1];
				x21 = x2[t];
				a2 = x2t1 + vt1;
				b2 = y2[t] + ut;

				z = z >= a ? z : a;
				z = z >= b ? z : b;
				z = z >= a2 ? z : a2;
				z = z >= b2 ? z : b2;
				z = z <= sc_mch_ ? z : sc_mch_;

				u_new = z - vt1;
				v_new = z - ut;
				a = a - (z - q_);
				b = b - (z - q_);
				a2 = a2 - (z - q2_);
				b2 = b2 - (z - q2_);

				if (a<0) {
					x_new = 0 - qe_;
				} else {
					x_new = a - qe_;
				}
				if (b<0) {
					y_new = 0 - qe_;
				} else {
					y_new = b - qe_;
				}
				if (a2<0) {
					x2_new = 0 - qe2_;
				} else {
					x2_new = a2 - qe2_;
				}
				if (b2<0) {
					y2_new = 0 - qe2_;
				} else {
					y2_new = b2 - qe2_;
				}

				if (t > st) {
					u[t-1] = u_next;
					v[t-1] = v_next;
					x[t-1] = x_next;
					y[t-1] = y_next;
					x2[t-1] = x2_next;
					y2[t-1] = y2_next;
				}
				if (t == en) {
					u[t] = u_new;
					v[t] = v_new;
					x[t] = x_new;
					y[t] = y_new;
					x2[t] = x2_new;
					y2[t] = y2_new;
					break;
				}
				u_next = u_new;
				v_next = v_new;
				x_next = x_new;
				y_next = y_new;
				x2_next = x2_new;
				y2_next = y2_new;
			}
		} 
		// !SCORE_ONLY && !RIGHT
		// !0x01 && !0x02
		else if (!(flag&KSW_EZ_RIGHT)) { // gap left-alignment
			uint8_t *pr = p + (size_t)r * n_col - st; // p is the matrix to store all the data when doing DP.
			off[r] = st, off_end[r] = en;
			int8_t u_new, v_new, x_new, y_new, x2_new, y2_new;
			int8_t u_next, v_next, x_next, y_next, x2_next, y2_next;
			for (t = st; t <= en; ++t) {
				int8_t d, z, a, b, a2, b2, xt1, x2t1, vt1, ut, tmp;

				z = s[t];
				xt1 = t == st ? x1_ : x[t-1];
				x1 = x[t];
				vt1 = t == st ? v1_ : v[t-1];
				v1 = v[t];
				a = xt1 + vt1;
				ut = u[t];
				b = y[t] + u[t];
				x2t1 = t == st ? x21_ : x2[t-1];
				x21 = x2[t];
				a2 = x2t1 + vt1;
				b2 = y2[t] + ut;

				if (a > z) {
					d = 1;
					z = a;
				} else {
					d = 0;
				}
				if (b > z) {
					d = 2;
					z = b;
				}
				if (a2 > z) {
					d = 3;
					z = a2;
				}
				if (b2 > z) {
					d = 4;
					z = b2;
				}
				if (sc_mch_ < z) {
					z = sc_mch_;
				}
				
				u_new = z - vt1;
				v_new = z - ut;
				a = a - (z - q_);
				b = b - (z - q_);
				a2 = a2 - (z - q2_);
				b2 = b2 - (z - q2_);

				if (a > 0) {
					x_new = a - qe_;
					d = d | 0x08;        // Set d[3] bit as high
				} else {
					x_new = 0 - qe_;
				}
				if (b > 0) {
					y_new = b - qe_;
					d = d | 0x10;        // Set d[4] as high
				} else {
					y_new = 0 - qe_;
				}
				if (a2 > 0) {
					x2_new = a2 - qe2_;
					d = d | 0x20;        // Set d[5] as high
				} else {
					x2_new = 0 - qe2_;
				}
				if (b2 > 0) {
					y2_new = b2 - qe2_;
					d = d | 0x40;        // Set d[6] as high
				} else {
					y2_new = 0 - qe2_;
				}
				pr[t] = (uint8_t)d;              // Store d

				if (t > st) {
					u[t-1] = u_next;
					v[t-1] = v_next;
					x[t-1] = x_next;
					y[t-1] = y_next;
					x2[t-1] = x2_next;
					y2[t-1] = y2_next;
				}
				if (t == en) {
					u[t] = u_new;
					v[t] = v_new;
					x[t] = x_new;
					y[t] = y_new;
					x2[t] = x2_new;
					y2[t] = y2_new;
					break;
				}
				u_next = u_new;
				v_next = v_new;
				x_next = x_new;
				y_next = y_new;
				x2_next = x2_new;
				y2_next = y2_new;

			}
		} 
		// !SCORE_ONLY && RIGHT
		// !0x01 && 0x02
		else { // gap right-alignment
			uint8_t *pr = p + (size_t)r * n_col - st; // p is the matrix to store all the data when doing DP.
			off[r] = st, off_end[r] = en; // off is a r_length*2 length int array, off_end is the pointer to the middle
			int8_t u_new, v_new, x_new, y_new, x2_new, y2_new;
			int8_t u_next, v_next, x_next, y_next, x2_next, y2_next;
			for (t = st; t <= en; ++t) {
				int8_t d, z, a, b, a2, b2, xt1, x2t1, vt1, ut, tmp;

				z = s[t];
				xt1 = t == st ? x1_ : x[t-1];
				x1 = x[t];
				vt1 = t == st ? v1_ : v[t-1];
				v1 = v[t];
				a = xt1 + vt1;
				ut = u[t];
				b = y[t] + u[t];
				x2t1 = t == st ? x21_ : x2[t-1];
				x21 = x2[t];
				a2 = x2t1 + vt1;
				b2 = y2[t] + ut;

				if (z>a)
				{
					d = 0;
					// z = z;
				}else
				{
					d = 1;
					z = a;
				}
				if (z<=b)
				{
					d = 2;
					z = b;
				}
				if (z<=a2)
				{
					d = 3;
					z = a2;
				}
				if (z<=b2)
				{
					d = 4;
					z = b2;
				}
				if (sc_mch_<z)
				{
					z = sc_mch_;
				}

				u_new = z - vt1;
				v_new = z - ut;
				a = a - (z - q_);
				b = b - (z - q_);
				a2 = a2 - (z - q2_);
				b2 = b2 - (z - q2_);

				if (a<0)
				{
					x_new = 0 - qe_;
				}else
				{
					x_new = a - qe_;
					d = d | 0x08;        // Set d[3] bit as high
				}
				if (b<0)
				{
					y_new = 0 - qe_;
				}else
				{
					y_new = b - qe_;
					d = d | 0x10;        // Set d[4] as high
				}
				if (a2<0)
				{
					x2_new = 0 - qe2_;
				}else
				{
					x2_new = a2 - qe2_;
					d = d | 0x20;        // Set d[5] as high
				}
				if (b2<0)
				{
					y2_new = 0 - qe2_;
				}else
				{
					y2_new = b2 - qe2_;
					d = d | 0x40;        // Set d[6] as high
				}
				pr[t] = (uint8_t)d;              // Store d
				

				if (t > st) {
					u[t-1] = u_next;
					v[t-1] = v_next;
					x[t-1] = x_next;
					y[t-1] = y_next;
					x2[t-1] = x2_next;
					y2[t-1] = y2_next;
				}
				if (t == en) {
					u[t] = u_new;
					v[t] = v_new;
					x[t] = x_new;
					y[t] = y_new;
					x2[t] = x2_new;
					y2[t] = y2_new;
					break;
				}
				u_next = u_new;
				v_next = v_new;
				x_next = x_new;
				y_next = y_new;
				x2_next = x2_new;
				y2_next = y2_new;

			} // Core loop ended (inner loop)
		}
		
		
		if (!approx_max) { // find the exact max with a 32-bit score array
			int32_t max_H, max_t;
			// compute H[], max_H and max_t
			if (r > 0) {
				int32_t HH[4], tt[4], en1 = en0, i;
				max_H = H[en0] = en0 > 0? H[en0-1] + u[en0] : H[en0] + v[en0]; // special casing the last element
				max_t = en0;
				for (t = st0; t < en0; t++) { // NOTE: why "-qe"? I don't think so.     this implements: H[t]+=v8[t]-qe; if(H[t]>max_H) max_H=H[t],max_t=t;			
					H[t]+=v8[t];
					if(H[t] > max_H){
						max_H = H[t];
						max_t = t;
					}
				}
			} else H[0] = v8[0] - qe, max_H = H[0], max_t = 0; // special casing r==0
			// update ez
			if (en0 == tlen - 1 && H[en0] > ez->mte)
				ez->mte = H[en0], ez->mte_q = r - en;
			if (r - st0 == qlen - 1 && H[st0] > ez->mqe)
				ez->mqe = H[st0], ez->mqe_t = st0;
			if (ksw_apply_zdrop(ez, 1, max_H, r, max_t, zdrop, e2)) {
				break;
			}
			if (r == qlen + tlen - 2 && en0 == tlen - 1)
				ez->score = H[tlen - 1];
		} else { // find approximate max; Z-drop might be inaccurate, too.
			if (r > 0) {
				if (last_H0_t >= st0 && last_H0_t <= en0 && last_H0_t + 1 >= st0 && last_H0_t + 1 <= en0) {
					int32_t d0 = v8[last_H0_t];
					int32_t d1 = u8[last_H0_t + 1];
					if (d0 > d1) H0 += d0;
					else H0 += d1, ++last_H0_t;
				} else if (last_H0_t >= st0 && last_H0_t <= en0) {
					H0 += v8[last_H0_t];
				} else {
					++last_H0_t, H0 += u8[last_H0_t];
				}
			} else H0 = v8[0] - qe, last_H0_t = 0;
			if ((flag & KSW_EZ_APPROX_DROP) && ksw_apply_zdrop(ez, 1, H0, r, last_H0_t, zdrop, e2)) {
				break;
			}
			if (r == qlen + tlen - 2 && en0 == tlen - 1)
				ez->score = H0;
		}
		
		last_st = st, last_en = en;
		
	}
	kfree(km, mem);
	if (!approx_max) kfree(km, H);
	if (with_cigar) { // backtrack
		int rev_cigar = !!(flag & KSW_EZ_REV_CIGAR);
		if (!ez->zdropped && !(flag&KSW_EZ_EXTZ_ONLY)) {
		ksw_backtrack(km, 1, rev_cigar, 0, (uint8_t*)p, off, off_end, n_col, tlen-1, qlen-1, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
		} else if (!ez->zdropped && (flag&KSW_EZ_EXTZ_ONLY) && ez->mqe + end_bonus > (int)ez->max) {
		ez->reach_end = 1;
		ksw_backtrack(km, 1, rev_cigar, 0, (uint8_t*)p, off, off_end, n_col, ez->mqe_t, qlen-1, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
		} else if (ez->max_t >= 0 && ez->max_q >= 0) {
		ksw_backtrack(km, 1, rev_cigar, 0, (uint8_t*)p, off, off_end, n_col, ez->max_t, ez->max_q, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
		}
		kfree(km, mem2); kfree(km, off);
	}
}
#endif // __SSE2__