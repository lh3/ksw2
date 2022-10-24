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

#ifdef DEBUG
#include <iostream>
FILE *align_score_file = NULL;
FILE *align_debug_file = NULL;
#endif

/* return non-zero if we won't see any mismatches */
__inline__ int ksw_check_par(const int8_t m, const int8_t *mat, const int8_t q,
              const int8_t e) {
    int max_sc, min_sc;
    for (int t = 1, max_sc = mat[0], min_sc = mat[1]; t < m * m; ++t) {
        max_sc = max_sc > mat[t] ? max_sc : mat[t];
        min_sc = min_sc < mat[t] ? min_sc : mat[t];
    }  // process the mat, get max_sc min_sc
    if (-min_sc > 2 * (q + e)) return 1;  // otherwise, we won't see any mismatches
    return 0;
}

/* t calculation:
    ref: even r :diag (i,j) = (0,0), (1,1), (2,2)... -> t = (n_col + 2)/2
         odd r  :diag j+1 (i,j+1) = (0,1), (1,2), (2,3)... -> t = (n_col + 2)/2
    anti-diag ref: i_ref = r/2, j_ref = (r+1)/2
*/
__inline__ int get_t(int i, int r, int n_col) {
    int offset = ((n_col + 2) / 2) - (r / 2);
    return i + offset;
}

__inline__ int get_i(int t, int r, int n_col) {
    int offset = ((n_col + 2) / 2) - (r / 2);
    return t - offset;
}

/* calculation z: max{ 	sc_{i,j},
						x_{i-1,j} + v_{i-1, j},
						y_{i,j-1}, u_{i,j-1},
						x2_{i-1,j} + v_{i-1,j},
						y2_{i,j-1} + u_{i,j-1}	}
	equation (5)
*/
__inline__ void ksw_cal_z(int8_t *z_ptr, int8_t sc, int8_t u, int8_t v, int8_t x,
                          int8_t y, int8_t x2, int8_t y2, const int8_t mm0){
    x = x + v;
    y = y + u;
    x2 = x2 + v;
    y2 = y2 + u;

    int8_t z = sc;
    z = x > z ? x : z;
    z = y > z ? y : z;
    z = x2 > z ? x2 : z;
    z = y2 > z ? y2 : z;
    z = z <= mm0 ? z : mm0;

    *z_ptr = z;
}

/* calculation z: max{ 	sc_{i,j},
                        x_{i-1,j} + v_{i-1, j},
                        y_{i,j-1}, u_{i,j-1},
                        x2_{i-1,j} + v_{i-1,j},
                        y2_{i,j-1} + u_{i,j-1}	}
    equation (5)
	return d for gap-left-alignment */
__inline__ int8_t ksw_cal_z_left_aligned(int8_t *z_ptr, int8_t sc, int8_t u, int8_t v, int8_t x, int8_t y,
			int8_t x2, int8_t y2, const int8_t mm0) {
    x = x + v;
    y = y + u;
    x2 = x2 + v;
    y2 = y2 + u;

	// find max
    int8_t z = sc;
    int8_t d = 0;

    if (x > z) {
        z = x;
        d = 1;
    }
	if (y > z) {
        z = y;
        d = 2;
    }
	if (x2 > z){
        z = x2;
        d = 3;
    }
	if (y2 > z){
        z = y2;
        d = 4;
    }
    z = z <= mm0 ? z : mm0;

    *z_ptr = z;

    return d;
}

/* calculation z: max{ 	sc_{i,j},
                        x_{i-1,j} + v_{i-1, j},
                        y_{i,j-1}, u_{i,j-1},
                        x2_{i-1,j} + v_{i-1,j},
                        y2_{i,j-1} + u_{i,j-1}	}
    equation (5)
	return d for gap-right-alignment */
__inline__ int8_t ksw_cal_z_right_aligned(int8_t *z_ptr, int8_t sc, int8_t u,
                                         int8_t v, int8_t x, int8_t y,
                                         int8_t x2, int8_t y2,
                                         const int8_t mm0) {
    x = x + v;
    y = y + u;
    x2 = x2 + v;
    y2 = y2 + u;

    // find max
    int8_t z = sc;
    int8_t d = 0;

    if (x >= z) {
        z = x;
        d = 1;
    }
    if (y >= z) {
        z = y;
        d = 2;
    }
    if (x2 >= z) {
        z = x2;
        d = 3;
    }
    if (y2 >= z) {
        z = y2;
        d = 4;
    }
    z = z <= mm0 ? z : mm0;

    *z_ptr = z;

    return d;
}

/* 	calculate x_ij, y_ij, x2_ij, y2_ij
	x_{i,j} = max{0, x_{i-1,j} + v_{i-1} - z_{ij} + q} - q - e}  equ(5)
*/
void ksw_cal_xy(int8_t *dst, int8_t x, const int8_t v, int8_t z,
                             const int8_t q, const int8_t e){
    x = x + v - z + q;
	if (x <= 0){
        x = 0;
    }
    *dst = x - q - e;
}

/* 	calculate x_ij, y_ij, x2_ij, y2_ij
    x_{i,j} = max{0, x_{i-1,j} + v_{i-1} - z_{ij} + q} - q - e}  equ(5)
    return 1 if a continuation on the E/F state (left aligned)
*/
int8_t ksw_cal_xy_left_aligned(int8_t *dst, int8_t x, const int8_t v, int8_t z,
                           const int8_t q, const int8_t e) {

    x = x + v - z + q;
    int d = 1;
    if (x <= 0) {
        x = 0;
        d = 0;
    }
    *dst = x - q - e;
    return d;
}

/* 	calculate x_ij, y_ij, x2_ij, y2_ij
    x_{i,j} = max{0, x_{i-1,j} + v_{i-1} - z_{ij} + q} - q - e}  equ(5)
    return 1 if a continuation on the E/F state (right aligned)
*/
int8_t ksw_cal_xy_right_aligned(int8_t *dst, int8_t x, const int8_t v,
                                        int8_t z, const int8_t q,
                                        const int8_t e) {
    x = x + v - z + q;
    int8_t d = 1;
    if (x < 0) {
        x = 0;
        d = 0;
    }
    *dst = x - q - e;
    return d;
}

template <int GEN_SC>
__inline__ int8_t ksw_cal_score( uint8_t target, uint8_t query, int8_t m, int8_t *mat, int8_t sc_N_){
    if (GEN_SC){
        if (query == m-1 || target == m-1) {
            return sc_N_;
        } else {
            return query == target ? mat[0] : mat[1];
        }
    } else {
        return mat[target * m + query];
    }
}


template <int SCORE_ONLY, int LEFT_ALIGNED=0>
int ksw_update_diag	(
	int8_t *sc, int8_t *u, int8_t *v, int8_t *x, int8_t *y, int8_t *x2, int8_t *y2,
    uint8_t *p, // indexed by t
    int32_t *H, int32_t *Hmax, int* rmax, // find max H
	const int r, const int t_st, const int t_en, int t_i_1, const int n_col, // start and end index in terms of t
	const int8_t q, const int8_t e, const int8_t q2, const int8_t e2, const int8_t mm0, const int8_t neta
){
    int8_t u_new=neta, v_new=neta, x_new=-q-e, y_new=-q-e, x2_new=-q2-e2, y2_new=-q2-e2;
    int32_t H_new = -q-e;
    for (int t = 0; t < n_col + 2; ++t, ++t_i_1) {
        // access previous matrix
        int8_t sc_elt, prev_u, prev_v, prev_x, prev_y, prev_x2, prev_y2;
        int32_t prev_H, prev_H_1;
        if (t >= t_st && t <= t_en){
            sc_elt = sc[t];           // s_{i,j}
            prev_u = u[t_i_1 + 1];    // u_{i,j+1}
            prev_v = v[t_i_1];        // v_{i-1,j}
            prev_x = x[t_i_1];        // x_{i-1,j}
            prev_y = y[t_i_1 + 1];    // y_{i,j-1}
            prev_x2 = x2[t_i_1];      // x2_{i-1,j}
            prev_y2 = y2[t_i_1 + 1];  // y2_{i,j-1}
            prev_H_1 = H[t_i_1];      // H_{i-1,j}
            prev_H = H[t_i_1 + 1];    // H_{i,j+1}
        }

        // update KZ matrix
        if (t >= 1){
            u[t-1] = u_new;
            v[t-1] = v_new;
            x[t-1] = x_new;
            y[t-1] = y_new;
            x2[t-1] = x2_new;
            y2[t-1] = y2_new;
            H[t-1] = H_new;
#ifdef DEBUG
                // if (t >= t_st && t <= t_en+1)
                //     fprintf(align_debug_file,
                //             "t %d, s %d, u %d, v %d, x %d,  y %d, x2 %d, y2
                //             %d H %d\n", t - 1, sc[t-1], u_new, v_new, x_new,
                //             y_new, x2_new, y2_new, H_new);
#endif
        }

        /* set default value */
        if (t < t_st || t > t_en) {
            u_new = v_new = neta;
            x_new = y_new = -q - e;
            x2_new = y2_new = -q2 - e2;
            H_new = KSW_NEG_INF;
            continue;
        }

        int8_t d, z;
        if (SCORE_ONLY){ // score only
            ksw_cal_z(&z, sc_elt, prev_u, prev_v, prev_x, prev_y, prev_x2, prev_y2, mm0);
        } else if (LEFT_ALIGNED){ // left aligned
            d = ksw_cal_z_left_aligned(&z, sc_elt, prev_u, prev_v, prev_x, prev_y,
                                       prev_x2, prev_y2, mm0);
        } else { // right aligned
            d = ksw_cal_z_right_aligned(&z, sc_elt, prev_u, prev_v, prev_x, prev_y,
                                        prev_x2, prev_y2, mm0);
        }
#ifdef DEBUG
        // fprintf(align_debug_file, "t %d t_{i-1} %d z %d prev_v %d prev_u %d\n", t, t_i_1, z, prev_v, prev_u);
#endif
        u_new = z - prev_v;
        v_new = z - prev_u;

        if (r == 0) {
            prev_H = -(int32_t)q-(int32_t)e;    // -q-e
        }

        if (SCORE_ONLY){
            ksw_cal_xy(&x_new, prev_x, prev_v, z, q, e);
            ksw_cal_xy(&y_new, prev_y, prev_u, z, q, e);
            ksw_cal_xy(&x2_new, prev_x2, prev_v, z, q2, e2);
            ksw_cal_xy(&y2_new, prev_y2, prev_u, z, q2, e2);
        } else if (LEFT_ALIGNED) {
            d = ksw_cal_xy_left_aligned(&x_new, prev_x, prev_v, z, q, e) ? d | 0x08 : d;
            d = ksw_cal_xy_left_aligned(&y_new, prev_y, prev_u, z, q, e) ? d | 0x10 : d;
            d = ksw_cal_xy_left_aligned(&x2_new, prev_x2, prev_v, z, q2, e2) ? d | 0x20 : d;
            d = ksw_cal_xy_left_aligned(&y2_new, prev_y2, prev_u, z, q2, e2) ? d | 0x40 : d;
            p[t] = (uint8_t)d;
        } else {
            d = ksw_cal_xy_right_aligned(&x_new, prev_x, prev_v, z, q, e) ? d | 0x08 : d;
            d = ksw_cal_xy_right_aligned(&y_new, prev_y, prev_u, z, q, e) ? d | 0x10 : d;
            d = ksw_cal_xy_right_aligned(&x2_new, prev_x2, prev_v, z, q2, e2) ? d | 0x20 : d;
            d = ksw_cal_xy_right_aligned(&y2_new, prev_y2, prev_u, z, q2, e2) ? d | 0x40 : d;
            p[t] = (uint8_t)d;
        }

        /* Calculate H anyways: for GPU */
        if (prev_H == KSW_NEG_INF) { // special casting the last element
            H_new = prev_H_1 + (int32_t)u_new;
        } else {
            H_new = prev_H + (int32_t)v_new;
        }

        // update H max & t max
        if (H_new > Hmax[t]){
            Hmax[t] = H_new;
            rmax[t] = r;
        }
#ifdef DEBUG
        // fprintf(align_debug_file, "%d|%d ", prev_H, H_new);
#endif
        }
    // update KZ matrix
    u[n_col+1] = u_new;
    v[n_col+1] = v_new;
    x[n_col+1] = x_new;
    y[n_col+1] = y_new;
    x2[n_col+1] = x2_new;
    y2[n_col+1] = y2_new;
    H[n_col+1] = H_new;
}

void ksw_extd2_cpp(
    /* Memory Space for kmalloc */
    void *km,
    /* Inputs */
    int qlen, const uint8_t *query, int tlen, const uint8_t *target,
    /* Parameters */
    int8_t m,
    const int8_t *mat,  // look up table. Constant across queries
    int8_t q, int8_t e, int8_t q2, int8_t e2,  // gap pennalty, constant
    int w,          // bandwidth, vary from query to query
    int zdrop,      // drop threshold, constant
    int end_bonus,  // ???, constant
    int flag,       // alignment flag, vary from query to query
    /* Output */
    ksw_extz_t *ez  // score and cigar
) {
    // debug
#ifdef DEBUG
    if (!align_debug_file) {
		align_debug_file = fopen("debug/test_sample_debug.output", "w+");
	}
    fprintf(align_debug_file, "H: \n");
    if (!align_score_file){
        align_score_file = fopen("debug/test_sample_score.output", "w+");
    }
    fprintf(align_score_file, "p: \n");
    // fprintf(align_debug_file, "q %d e %d q2 %d e2 %d\n", q, e, q2, e2);
    // fprintf(align_debug_file,
    //         "(r, t, i | u, v, x, y, x2, y2, H, Hmax, rmax, p)\n");
#endif

    /* Score Generation parameters */
    int long_thres, long_diff; // derived from q, e.
    int8_t sc_N_; // derived from mat and e2

    int min_n_col;  // min of qlen, tlen, w, number of col to record.
    int n_col;  // min { w, max{qlen, tlen}} 

    /* Score generation intermediate vals */
    int8_t *u, *v, *x, *y, *x2, *y2; // Suzuki-Kasahara formulation
	int8_t *sc;	// score
    int32_t *H, *Hmax;  // intermediate val for finding max after score generation
    int *rmax; // val for finding max location after score generation

    /* options */
    int with_cigar = !(flag & KSW_EZ_SCORE_ONLY);
    assert( !(flag & KSW_EZ_APPROX_MAX) ); // GPU doesn't support max_approx

    /* Score Generation output for backtracking */
    int *off = 0,
        *off_end = 0;  // points to the start and end of
    uint8_t *p = 0;

    ksw_reset_extz(ez);
	if (m <= 1 || qlen <= 0 || tlen <= 0) return;

    int tmp;
    if (q2 + e2 < q + e) tmp = q, q = q2, q2 = tmp, tmp = e, e = e2, e2 = tmp; // make sure q+e no larger than q2+e2

    if (ksw_check_par(m, mat, q, e)) return;

    /* Derive Score intermediate values */
	//long_thres
    long_thres = e != e2? (q2 - q) / (e - e2) - 1 : 0;
	if (q2 + e2 + long_thres * e2 > q + e + long_thres * e)
		++long_thres;
	//long_diff
	long_diff = long_thres * (e - e2) - (q2 - q) - e2;
    // sc_N_
    sc_N_ = mat[m * m - 1] == 0 ? -e2 : mat[m * m - 1];

    /* memory size  */
    if (w < 0) w = tlen > qlen ? tlen : qlen;
    n_col = qlen > tlen ? qlen : tlen;
    n_col = n_col < (w + 1) ? n_col : w + 1;

    min_n_col = qlen < tlen ? qlen : tlen;
    min_n_col = min_n_col < (w + 1) ? min_n_col : w + 1;

    // int full_target = 0;
    // if (n_col == qlen) {
    //     full_query = 1;
    //     n_col = tlen;
    // }
    // int real_n_col = qlen < n_col ? qlen : n_col;

    /* Allocate memory & initialize intermediate value */
	// n_col + 2 is enough for u, v, x, y, x2, y2, s
    /* u, v, x, y, x2, y2, sc, H, Hmax, rmax INDEXING (by t):
        t = 0: Reserved for initialization value. Always outside of the band!
        t = 1 - w-1/w: within band.
        t = w: maybe within band or initialization value.

        Details about t to i tranlation: see get_i & get_t
    */
    u = (int8_t *)kmalloc(km, n_col + 2);
    memset(u, -q - e, n_col+2);
    v = (int8_t *)kmalloc(km, n_col + 2);
    memset(v, -q - e, n_col+2);
    x = (int8_t *)kmalloc(km, n_col + 2);
    memset(x, -q - e, n_col+2);
    y = (int8_t *)kmalloc(km, n_col + 2);
    memset(y, -q - e, n_col+2);
    x2 = (int8_t *)kmalloc(km, n_col + 2);
    memset(x2, -q2 - e2, n_col+2);
    y2 = (int8_t *)kmalloc(km, n_col + 2);
    memset(y2, -q2 - e2, n_col+2);
    sc = (int8_t *)kmalloc(km, n_col + 2);

    H = (int32_t*)kmalloc(km, (n_col+2) * sizeof(int32_t));
    for (int t = 0; t <= n_col; ++t) H[t] = KSW_NEG_INF;
    Hmax = (int32_t *)kmalloc(km, (n_col + 2) * sizeof(int32_t));
    for (int t = 0; t <= n_col; ++t) Hmax[t] = KSW_NEG_INF; // Hmax[i] = -inf
    rmax = (int *)kmalloc(km, (n_col + 2) * sizeof(int));
    for (int t = 0; t <= n_col; ++t) rmax[t] = -1;

    if (with_cigar) {
        size_t p_size = ((size_t)qlen + tlen - 1) * (size_t)min_n_col + 1;
        p = (uint8_t *)kmalloc(km, p_size);  // qlen + tlen - 1 = number of diagonal
        assert(p != NULL);
        // In the backtrack matrix, value p[] has the following structure:
		//   bit 0-2: which type gets the max - 0 for H, 1 for E, 2 for F, 3 for \tilde{E} and 4 for \tilde{F}
		//   bit 3/0x08: 1 if a continuation on the E state (bit 5/0x20 for a continuation on \tilde{E})
		//   bit 4/0x10: 1 if a continuation on the F state (bit 6/0x40 for a continuation on \tilde{F})
		off = (int*)kmalloc(km, (qlen + tlen - 1) * sizeof(int) * 2); // off is a int array using r to index
		off_end = off + qlen + tlen - 1;
	}

#ifdef DEBUG
    if (!with_cigar)
        fprintf(align_debug_file, "non score\n");
    else if (flag & KSW_EZ_RIGHT)
        fprintf(align_debug_file, "right aligned\n");
    else
        fprintf(align_debug_file, "left aligned\n");
#endif

    int t_st, t_en;
    for (int r = 0; r < qlen + tlen - 1; ++r) {   // r: iterate through anti-diag

        /* NOTE: find anti-diag boundaries (in terms of i) */
        int st = 0, en = tlen-1;
        if (st < r - qlen + 1) st = r - qlen + 1;  // Correct the value of st iif r is large, during each iteration (Can be explained using figure)
		if (en > r) en = r; // Correct the en if r is small
		if (st < (r-w+1)>>1) st = (r-w+1)>>1; // take the ceil, choose the band position
		if (en > (r+w)>>1) en = (r+w)>>1; // take the floor, choose the band position
		//DEBUG: change to assert
        #ifdef DEBUG
		if (st > en) {
			ez->zdropped = 1;
			printf("break due to st > en \n"); // debug
			break;
		}
        #endif

        // printf("r %d st %d en %d n_col %d t_st %d t_en %d\n", r, st, en, n_col, t_st, t_en);
        t_st = get_t(st, r, n_col);
        t_en = get_t(en, r, n_col);
        assert(t_st > 0);

        int t_i_1 = - (r & 0x1);
        assert(t_i_1 + t_st >= 0);
        assert(t_i_1 + t_en + 1 < n_col + 2);

        /* NOTE: find boundary values for neta */
        int8_t neta = r + 1 < long_thres    ? -e
                    : r + 1 == long_thres   ? long_diff
                                            : -e2;

		/* NOTE: Preprocess score */
		// match/mismatch only
		if (!(flag & KSW_EZ_GENERIC_SC)) {
			for (int i = st, t = t_st; i <= en; t++, i++) {
				int8_t query_elt, target_elt, score;
                query_elt = target[i];  		// base on reference chain
                target_elt = query[r - i];  	// base on query chain
                if (query_elt == m - 1 || target_elt == m - 1) {
                    score = sc_N_;
                } else {
                    score = query_elt == target_elt ? mat[0] : mat[1];
                }
                *(sc+t) = score;
            }
		} else { // score
            for (int i = st, t = t_st; i <= en; t++, i++){
                ((uint8_t*)sc)[t] = mat[target[i] * m + query[r-i]];
				// Preprocess the score, matrix is like a lookup table.
            }
		}
        
        if (with_cigar){
            off[r] = st;
            off_end[r] = en;
        }
        // DEBUG: check p access
        uint8_t *pr = p + r * min_n_col - t_st;
        assert(pr + t_en < p + (r+1) * min_n_col);
        assert(pr + t_st == p + r * min_n_col);

        /* NOTE: update ksw matrix & H, p */
        if (!with_cigar) {  // score only
            ksw_update_diag<1>(sc, u, v, x, y, x2, y2,
                        pr, H,
                        Hmax, rmax,
	    				r, t_st, t_en, t_i_1, n_col,
	    				q, e, q2, e2, mat[0], neta);
        } else if (!(flag & KSW_EZ_RIGHT)) {  // gap left_aligned
            ksw_update_diag<0, 1>(sc, u, v, x, y, x2, y2,
                        pr, H,
                        Hmax, rmax,
	    				r, t_st, t_en, t_i_1, n_col,
	    				q, e, q2, e2, mat[0], neta);
        } else {  // right algined
            ksw_update_diag<0, 0>(sc, u, v, x, y, x2, y2,
                        pr, H,
                        Hmax, rmax,
	    				r, t_st, t_en, t_i_1, n_col,
	    				q, e, q2, e2, mat[0], neta);
	    }

        /* NOTE: update mte & mqe */
        if (en == tlen - 1 && H[t_en] > ez->mte){
            ez->mte = H[t_en], ez->mte_q = r - en;
        }
        if (r - st == qlen - 1 && H[t_st] > ez->mqe){
            ez->mqe = H[t_st], ez->mqe_t = st;
        }

        // DEBUG: debug output
#ifdef DEBUG
        fprintf(align_debug_file, "#%d (st=%d en=%d) ", r, st, en);
        for (int t = t_st; t <= t_en; ++t) {
        //     if (!align_debug_file) {
        //         align_debug_file = fopen("debug/test_sample_debug.output", "w+");
        //     }
            fprintf(align_debug_file, "%d ", H[t]);
        //     fprintf(align_debug_file, "(%d,%d,%d|%d,%d,%d,%d,%d,%d,%d,%d,%d,0x%x)\n",
        //             r, t, get_i(t, r, n_col), ((int8_t *)u)[t], ((int8_t *)v)[t],
        //             ((int8_t *)x)[t], ((int8_t *)y)[t], ((int8_t *)x2)[t],
        //             ((int8_t *)y2)[t], ((int32_t *)H)[t], ((int32_t *)Hmax)[t],
        //             ((int *)rmax)[t], p[r*min_n_col - t_st + t]);  // for debugging
        //     if (!align_score_file) {
        //         align_score_file = fopen("debug/test_sample_score.output", "w+");
        //         fprintf(align_score_file, "(r, t | u, v, x, y)\n");
        //     }
        //     fprintf(align_score_file, "(%d,%d|%d,%d,%d,%d)", 
        //         r, t-t_st+st, ((int8_t*)u)[t], ((int8_t*)v)[t], ((int8_t*)x)[t], 
        //         ((int8_t*)y)[t]); // for debugging
        }
        // fprintf(align_score_file, "\n");
        fprintf(align_debug_file, "\n");
#endif
    }  // r loop NOTE: output of the loop: Hmax, rmax, ez, p

    // NOTE: find max for ez. When multiple entries has the same score, always chose the one with smaller r. 
    int max_r = -1;
    for (int t = 1; t <= n_col; ++t){
        if (Hmax[t] > (int32_t)ez->max || (Hmax[t] == (int32_t)ez->max && rmax[t] < max_r)) {
            max_r = rmax[t];
            ez->max = Hmax[t];
            ez->max_t = get_i(t, rmax[t], n_col);   // max_t
            ez->max_q = rmax[t] - ez->max_t;        // r - max_t
        }
    }

    ez->score = H[t_en];

#ifdef DEBUG

    // fprintf(
    //     align_debug_file,
    //     "ez: max %d zdropped %d max_q %d max_t %d mte %d mte_q %d score %d\n",
    //     ez->max, ez->zdropped, ez->max_q, ez->max_t, ez->mte, ez->mte_q,
    //     ez->score);
    for (int i = 0; i < qlen + tlen - 1; i++) {
        int len = i + 1 < qlen + tlen - i - 1 ? i + 1 : qlen + tlen - i - 1;
        len = min_n_col < len ? min_n_col : len;
        fprintf(align_score_file, "#%d ", i);
        for (int j = 0; j < len; j++)
            fprintf(align_score_file, "%x ", ((uint8_t *)p)[i * min_n_col + j]);
        fprintf(align_score_file, "\n");
    }

    // for (int i = 0; i < qlen + tlen - 1; i++) {
    //     fprintf(align_debug_file, "%d ", off[i]);
    // }fprintf(align_debug_file, "\n");
    // for (int i = 0; i < qlen + tlen - 1; i++) {
    //     fprintf(align_debug_file, "%d ", off_end[i]);
    // }fprintf(align_debug_file, "\n");

#endif

    kfree(km, u);
    kfree(km, v);
    kfree(km, x);
    kfree(km, y);
    kfree(km, x2);
    kfree(km, y2);
    kfree(km, sc);
    kfree(km, H);
    kfree(km, Hmax);
    kfree(km, rmax);
    if (with_cigar) { // backtrack
                      //
		int rev_cigar = !!(flag & KSW_EZ_REV_CIGAR);
		if (!ez->zdropped && !(flag&KSW_EZ_EXTZ_ONLY)) {
		ksw_backtrack(km, 1, rev_cigar, 0, (uint8_t*)p, off, off_end, min_n_col, tlen-1, qlen-1, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
		} else if (!ez->zdropped && (flag&KSW_EZ_EXTZ_ONLY) && ez->mqe + end_bonus > (int)ez->max) {
		ez->reach_end = 1;
		ksw_backtrack(km, 1, rev_cigar, 0, (uint8_t*)p, off, off_end, min_n_col, ez->mqe_t, qlen-1, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
		} else if (ez->max_t >= 0 && ez->max_q >= 0) {
		ksw_backtrack(km, 1, rev_cigar, 0, (uint8_t*)p, off, off_end, min_n_col, ez->max_t, ez->max_q, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
		}
		kfree(km, p); kfree(km, off);
	}
}
#endif // __SSE2__
