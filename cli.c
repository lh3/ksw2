#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "ksw2.h"

unsigned char seq_nt4_table[256] = {
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

static void ksw_gen_simple_mat(int m, int8_t *mat, int8_t a, int8_t b)
{
	int i, j;
	a = a < 0? -a : a;
	b = b > 0? -b : b;
	for (i = 0; i < m - 1; ++i) {
		for (j = 0; j < m - 1; ++j)
			mat[i * m + j] = i == j? a : b;
		mat[i * m + m - 1] = 0;
	}
	for (j = 0; j < m; ++j)
		mat[(m - 1) * m + j] = 0;
}

int main(int argc, char *argv[])
{
	int8_t a = 1, b = 1, q = 1, e = 1;
	int c, i, qlen, tlen, n_cigar, score;
	char *qseq, *tseq, *algo = "gg2";
	uint32_t *cigar;
	int8_t mat[25];
	ksw_extz_t ez;

	while ((c = getopt(argc, argv, "t:")) >= 0) {
		if (c == 't') algo = optarg;
	}
	if (argc - optind < 2) {
		fprintf(stderr, "Usage: ksw2-global <DNA-target> <DNA-query>\n");
		return 1;
	}
	memset(&ez, 0, sizeof(ksw_extz_t));
	ksw_gen_simple_mat(5, mat, a, -b);
	tseq = argv[optind], qseq = argv[optind+1];
	tlen = strlen(tseq);
	qlen = strlen(qseq);
	for (i = 0; i < qlen; ++i)
		qseq[i] = seq_nt4_table[(uint8_t)qseq[i]];
	for (i = 0; i < tlen; ++i)
		tseq[i] = seq_nt4_table[(uint8_t)tseq[i]];
	if (strcmp(algo, "gg") == 0)           score = ksw_gg(0, qlen, (uint8_t*)qseq, tlen, (uint8_t*)tseq, 5, mat, q, e, qlen > tlen? qlen : tlen, &n_cigar, &cigar);
	else if (strcmp(algo, "gg2") == 0)     score = ksw_gg2(0, qlen, (uint8_t*)qseq, tlen, (uint8_t*)tseq, 5, mat, q, e, qlen > tlen? qlen : tlen, &n_cigar, &cigar);
	else if (strcmp(algo, "gg2_sse") == 0) score = ksw_gg2_sse(0, qlen, (uint8_t*)qseq, tlen, (uint8_t*)tseq, 5, mat, q, e, qlen > tlen? qlen : tlen, &n_cigar, &cigar);
	else if (strcmp(algo, "extz_sse") == 0) {
		ksw_extz_sse(0, qlen, (uint8_t*)qseq, tlen, (uint8_t*)tseq, 5, mat, q, e, qlen > tlen? qlen : tlen, 100, 0, &ez);
		n_cigar = ez.n_cigar, cigar = ez.cigar, score = ez.score;
		printf("max: %d; (%d,%d)\n", ez.max, ez.max_t, ez.max_q);
		printf("mqe: %d; %d\n", ez.mqe, ez.mqe_t);
	} else abort();
	printf("%d\t", score);
	for (i = 0; i < n_cigar; ++i)
		printf("%d%c", cigar[i]>>4, "MID"[cigar[i]&0xf]);
	printf("\n");
	return 0;
}
