#ifndef KSW2_H_
#define KSW2_H_

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Global alignment with moving a band
 *
 * @param km        memory pool, when used with kalloc
 * @param qlen      query length
 * @param query     query sequence with 0 <= query[i] < m
 * @param tlen      target length
 * @param target    target sequence with 0 <= target[i] < m
 * @param m         number of residue types
 * @param mat       m*m scoring mattrix in one-dimension array
 * @param gapo      gap open penalty; a gap of length l cost "-(gapo+l*gape)"
 * @param gape      gap extension penalty
 * @param w         band width
 * @param n_cigar   (out) number of CIGAR elements
 * @param cigar     (out) BAM-encoded CIGAR; caller need to deallocate with kfree(km, )
 *
 * @return          score of the alignment
 */
int ksw_gg(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t gapo, int8_t gape, int w, int *n_cigar_, uint32_t **cigar_);
int ksw_gg2(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t gapo, int8_t gape, int w, int *n_cigar_, uint32_t **cigar_);
int ksw_gg2_sse(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t gapo, int8_t gape, int w, int *n_cigar_, uint32_t **cigar_);

#ifdef __cplusplus
}
#endif

/********************************
 * Private macros and functions *
 ********************************/

#ifdef HAVE_KALLOC
#include "kalloc.h"
#else
#include <stdlib.h>
#define kmalloc(km, size) malloc((size))
#define kcalloc(km, count, size) calloc((count), (size))
#define krealloc(km, ptr, size) realloc((ptr), (size))
#define kfree(km, ptr) free((ptr))
#endif

static inline uint32_t *ksw_push_cigar(void *km, int *n_cigar, int *m_cigar, uint32_t *cigar, int op, int len)
{
	if (*n_cigar == 0 || op != (cigar[(*n_cigar) - 1]&0xf)) {
		if (*n_cigar == *m_cigar) {
			*m_cigar = *m_cigar? (*m_cigar)<<1 : 4;
			cigar = krealloc(km, cigar, (*m_cigar) << 2);
		}
		cigar[(*n_cigar)++] = len<<4 | op;
	} else cigar[(*n_cigar)-1] += len<<4;
	return cigar;
}

#endif
