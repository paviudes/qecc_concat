#ifndef DECODE_H
#define DECODE_H

// Compute the probabilities of logical classes, given the probabilities of all Pauli errors.
extern void ComputeCosetProbsLevelOne(double *pauli_probs, int nlogs, int nstabs, double **cosetprobs);

// Compute the probabilities of logical classes, given the probabilities of logical I, X, Y and Z of the previous level.
extern void ComputeCosetProbs(int synd, double **pauli_probs, int ****TLS, int nphys, int nlogs, int nstabs, double *cosetprobs);

// Compute the logical whose coset probability is the highest.
extern int ArgMax(double *cosetprobs, int nlogs);

#endif /* DECODE_H */