#ifndef DECODE_H
#define DECODE_H

// Compute the probabilities of logical classes, given the probabilities of all Pauli errors.
extern void ComputeCosetProbsLevelOne(double *pauli_probs, int nlogs, int nstabs, double **cosetprobs);

// Compute the probabilities of logical classes, given the probabilities of logical I, X, Y and Z of the previous level.
extern void ComputeCosetProbs(int synd, double **pauli_probs, int ****LST, int nphys, int nlogs, int nstabs, double *cosetprobs);

// Compute the logical whose coset probability is the highest.
extern int ArgMax(double *cosetprobs, int nlogs);

/*
	Rotate a probability of Pauli operators by a Pauli.
	Letting I --> 0, X --> 1, Y --> 2 and Z --> 3
	I: [0, 1, 2, 3]
	X: [1, 0, 3, 2]
	Y: [2, 3, 0, 1]
	Z: [3, 2, 1, 0]
*/
extern void RotatePauli(double *arr, int size, int pauli);

#endif /* DECODE_H */