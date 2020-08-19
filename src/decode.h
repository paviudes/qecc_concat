#ifndef DECODE_H
#define DECODE_H

// Compute the probabilities of logical classes, given the probabilities of I, X, Y and Z.
extern void ComputeCosetProbs(double **pauli_probs, int ****TLS, int nphys, int nlogs, int nstabs, double **cosetprobs);

// Compute the logical whose coset probability is the highest.
extern void GetMaxCoset(double **cosetprobs, int nstabs, int nlogs, int *corrections);

#endif /* DECODE_H */