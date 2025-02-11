#ifndef QECC_H
#define QECC_H

#include <complex.h>
#include "memory.h"

struct qecc_t {
  int N;
  int K;
  int D;
  int nlogs;
  int nstabs;
  int **restrict projector;
  int ***restrict action;
  double complex **restrict phases;
  int *dclookup; // Lookup table specifying the corrections for each syndrome.
  int ****LST;
};

// Allocate memory for the elements of the quantum error correcting code.
extern void InitQECC(struct qecc_t *qecc);

// Allocate memory allocated to the elements of the quantum error correcting
// code..
extern void FreeQECC(struct qecc_t *qecc);

// Explicitly set the full process matrix of a correlated channel
extern void SetFullProcessMatrix(struct qecc_t *qecc, struct simul_t *sim, double *process, int isPauli);

// Compute the effective logical channel, when error correction is applied over
// a set of input physical channels.
extern void SingleShotErrorCorrection(int isPauli, int iscorr, int dcalg, int frame, struct qecc_t *qecc, struct simul_t *sim, struct constants_t *consts, int is_cosetprobs_computed);

#endif /* QECC_H */
