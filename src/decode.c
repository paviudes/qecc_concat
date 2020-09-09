#include <stdlib.h>
#include <stdio.h>
#include "printfuns.h"
#include "linalg.h"
#include "decode.h"
#include <math.h>

/*
	1. Remove coset_sum from ComputeCosetProbs and ComputeCosetProbsLevelOne functions.
*/

void ComputeCosetProbs(int synd, double **pauli_probs, int ****LST, int nphys, int nlogs, int nstabs, double *cosetprobs){
	// Compute the probabilty of the 4^k logical classes for each syndrome.
	// PrintDoubleArray2D(pauli_probs, "Pauli probs to the coset probs computation.", nphys, nlogs);
	int q, s, l;
	double prob;
	// double cosetprobs_sum = 0;
	for (l = 0; l < nlogs; l ++){
		cosetprobs[l] = 0;
		for (s = 0; s < nstabs; s ++){
			prob = 1;
			for (q = 0; q < nphys; q ++)
				prob *= pauli_probs[q][LST[l][s][synd][q]];
			cosetprobs[l] += prob;
		}
		// cosetprobs_sum += cosetprobs[synd][l];
	}
	Normalize(cosetprobs, nlogs);
	// printf("Coset probabilities for s = %d.\n", synd);
	// PrintDoubleArray1D(cosetprobs, "Coset probabilities", nlogs);
}

void ComputeCosetProbsLevelOne(double *pauli_probs, int nlogs, int nstabs, double **cosetprobs){
	// Compute the probabilty of the 4^k logical classes for each syndrome, for level 1.
	// PrintDoubleArray1D(pauli_probs, "Pauli probs for the level 1 coset probs computation.", nlogs * nstabs * nstabs);
	int synd, s, l;
	// double cosetprobs_sum = 0;
	for (synd = 0; synd < nstabs; synd ++){
		for (l = 0; l < nlogs; l ++){
			cosetprobs[synd][l] = 0;
			for (s = 0; s < nstabs; s ++)
				cosetprobs[synd][l] += pauli_probs[l * nstabs * nstabs + s * nstabs + synd];
			// cosetprobs_sum += cosetprobs[synd][l];
		}
		// PrintDoubleArray1D(cosetprobs[synd], "Before normalization", nlogs);
		Normalize(cosetprobs[synd], nlogs);
		// PrintDoubleArray1D(cosetprobs[synd], "After normalization", nlogs);
		// printf("Level one coset probabilities for s = %d.\n", synd);
		// PrintDoubleArray1D(cosetprobs[synd], "Coset probabilities", nlogs);
	}
	// printf("Sum of level one coset probabilities for all syndromes = %g.\n", cosetprobs_sum);
}

int ArgMax(double *arr, int size){
	// Compute the coset with maximum probability, for each syndrome.
	double prec = pow(10,12);
	int i, amax = 0;
	for (i = 1; i < size; i ++)
		if ( floorf(arr[i] * prec) / prec > floorf(arr[amax] * prec) / prec)
			amax = i;
	return amax;
}
