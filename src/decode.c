#include <stdlib.h>
#include <stdio.h>
#include "printfuns.h"
#include "decode.h"

void ComputeCosetProbs(double **pauli_probs, int ****TLS, int nphys, int nlogs, int nstabs, double **cosetprobs){
	// Compute the probabilty of the 4^k logical classes for each syndrome.
	// PrintDoubleArray2D(pauli_probs, "Pauli probs to the coset probs computation.", nphys, nlogs);
	int synd, q, s, l;
	double prob, cosetprobs_sum = 0;
	for (synd = 0; synd < nstabs; synd ++){
		for (l = 0; l < nlogs; l ++){
			cosetprobs[synd][l] = 0;
			for (s = 0; s < nstabs; s ++){
				prob = 1;
				for (q = 0; q < nphys; q ++)
					prob *= pauli_probs[q][TLS[synd][l][s][q]];
				cosetprobs[synd][l] += prob;
			}
			cosetprobs_sum += cosetprobs[synd][l];
		}
		// printf("Coset probabilities for s = %d.\n", synd);
		// PrintDoubleArray1D(cosetprobs[synd], "Coset probabilities", nlogs);
	}
	// printf("Sum of coset probabilities for all syndromes = %g.\n", cosetprobs_sum);
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