#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>
#include "ksw2.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

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

static void global_aln(const char *algo, void *km, char *qseq, char *tseq, int8_t m, const int8_t *mat, int8_t q, int8_t e, int w, int zdrop, int flag, ksw_extz_t *ez)
{
	int i, qlen, tlen;
	qlen = strlen(qseq);
	tlen = strlen(tseq);
	for (i = 0; i < qlen; ++i)
		qseq[i] = seq_nt4_table[(uint8_t)qseq[i]];
	for (i = 0; i < tlen; ++i)
		tseq[i] = seq_nt4_table[(uint8_t)tseq[i]];
	if (w < 0) w = qlen > tlen? qlen : tlen;
	if (strcmp(algo, "gg") == 0)           ez->score = ksw_gg(0, qlen, (uint8_t*)qseq, tlen, (uint8_t*)tseq, 5, mat, q, e, w, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
	else if (strcmp(algo, "gg2") == 0)     ez->score = ksw_gg2(0, qlen, (uint8_t*)qseq, tlen, (uint8_t*)tseq, 5, mat, q, e, w, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
	else if (strcmp(algo, "gg2_sse") == 0) ez->score = ksw_gg2_sse(0, qlen, (uint8_t*)qseq, tlen, (uint8_t*)tseq, 5, mat, q, e, w, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
	else if (strcmp(algo, "extz2_sse") == 0) ksw_extz2_sse(0, qlen, (uint8_t*)qseq, tlen, (uint8_t*)tseq, 5, mat, q, e, w, zdrop, flag, ez);
	else abort();
}

int main(int argc, char *argv[])
{
	int8_t a = 1, b = 1, q = 1, e = 1;
	int c, i, pair = 1, w = -1;
	char *algo = "gg2";
	int8_t mat[25];
	ksw_extz_t ez;
	gzFile fp[2];

	while ((c = getopt(argc, argv, "t:w:")) >= 0) {
		if (c == 't') algo = optarg;
		else if (c == 'w') w = atoi(optarg);
	}
	if (argc - optind < 2) {
		fprintf(stderr, "Usage: ksw2-global <DNA-target> <DNA-query>\n");
		return 1;
	}
	memset(&ez, 0, sizeof(ksw_extz_t));
	ksw_gen_simple_mat(5, mat, a, -b);
	fp[0] = gzopen(argv[optind], "r");
	fp[1] = gzopen(argv[optind+1], "r");

	if (fp[0] == 0 && fp[1] == 0) {
		global_aln(algo, 0, argv[optind+1], argv[optind], 5, mat, q, e, w, 100, 0, &ez);
		printf("%d\t", ez.score);
		for (i = 0; i < ez.n_cigar; ++i)
			printf("%d%c", ez.cigar[i]>>4, "MID"[ez.cigar[i]&0xf]);
		printf("\n");
	} else if (fp[0] && fp[1]) {
		kseq_t *ks[2];
		ks[0] = kseq_init(fp[0]);
		ks[1] = kseq_init(fp[1]);
		if (pair) {
			while (kseq_read(ks[0]) > 0) {
				if (kseq_read(ks[1]) <= 0) break;
				global_aln(algo, 0, ks[0]->seq.s, ks[1]->seq.s, 5, mat, q, e, w, 100, 0, &ez);
				printf("%s\t%s\t%d\t", ks[0]->name.s, ks[1]->name.s, ez.score);
				for (i = 0; i < ez.n_cigar; ++i)
					printf("%d%c", ez.cigar[i]>>4, "MID"[ez.cigar[i]&0xf]);
				printf("\n");
			}
		}
		kseq_destroy(ks[0]);
		kseq_destroy(ks[1]);
	}
	return 0;
}
