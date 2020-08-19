#include <stdlib.h>
#include "qcode.h"
#include "decode.h"

void ComputeCosetProbs(struct qecc_t qcode, double **pauli_probs, double **cosetprobs){
	// Compute the probabilty of the 4^k logical classes for each syndrome.
	int synd, s, l;
	for (synd = 0; synd < qcode->nlogs; synd ++){
		for (l = 0; l < qcode->nlogs; l ++){
			cosetprobs[synd][l] = 0;
			for (s = 0; s < qcode->nstabs; synd ++){
				prob = 1;
				for (q = 0; q < qcode->N; q ++)
					prob *= pauli_probs[q][(qcode->TLS)[synd][l][s][q]];
				cosetprobs[synd][l] += prob;
			}
		}
	}
}

void GetMaxCoset(double **cosetprobs, int nstabs, int nlogs, int *corrections){
	// Compute the coset with maximum probability, for each syndrome.
	int l, s;
	double max = 0;
	for (s = 0; s < nstabs; s ++){
		max = 0;
		for (l = 0; l < nlogs; l ++){
			if (cosetprobs[s][l] > max){
				max = cosetprobs[s][l];
				corrections[s] = l;
			}
		}
	}
}