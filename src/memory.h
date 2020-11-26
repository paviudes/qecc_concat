#ifndef MEMORY_H
#define MEMORY_H

#include <complex.h>

struct simul_t
{
  // Physical noise processes.
  char *restrict chname;
  int iscorr;
  double *restrict physical;
  double *restrict mpinfo;
  double ***restrict virtchan;
  double **restrict pauli_probs;
  // Quantum error correction and logical noise
  double **restrict tpnorms;
  double ***restrict logical;
  double ****restrict process;
  double *restrict syndprobs;
  double *restrict cumulative;
  double *restrict levelOneSynds;
  double *restrict levelOneCumul;
  int *restrict corrections;
  double ***restrict effprocess;
  double **restrict cosetprobs;
  double ***restrict levelOneChannels;
  double **restrict levelOneCosets;
  int *restrict frames;
  // Metrics.
  int nmetrics;
  char **restrict metricsToCompute;
  double **restrict metricValues;
  // Syndrome sampling.
  long nstats;
  int nlevels;
  int maxbin;
  int importance;
  double *restrict levelOneImpDist;
  double *restrict levelOneImpCumul;
  double infidelity;
  double *restrict outlierprobs;
  double **restrict sampling;
  long *restrict statsperlevel;
  int nbins;
  int ****restrict bins;
  double **restrict sumsq;
  double **restrict variance;
  int nbreaks;
  long *restrict runstats;
  double **restrict runavg;
  double threshold;
  // Decoder
  int *decoders;
  int hybrid;
  int **restrict decbins;
  int *restrict ndecbins;
  // Randomized compiling of QEC gates
  int rc;
  int *rcpauli;
};

// Initialize the elements that pertain to the montecarlo simulation of
// channels.
extern void AllocSimParams(struct simul_t *simul, int nphys, int nenc);

// Free memory allocated to the elements of the simulation structure that do not
// depend on the QECC.
extern void FreeSimParams(struct simul_t *simul, int nphys, int nenc);

/*
        Allocate memory for parameters of the simulation structure that depend
   upon the QECC used. This memeory will be reallocated everytime there is a new
   QECC, i.e, at a new concatenation level. These parameters are virtual,
   logical, syndprobs, cumulative, levelZeroSynds, levelZeroCumulative,
   levelZeroImpDist, levelZeroImpCumul, process, corrections, effprocess,
   effective, levelZeroChannels
*/
extern void AllocSimParamsQECC(struct simul_t *simul, int nphys, int nenc);

// Free the memory allocated to simulation parameters that depend on the QECC.
extern void FreeSimParamsQECC(struct simul_t *simul, int nphys, int nenc);

// Determine the number of independent logical channels at every level, that
// determine the logical channels of higher levels.
extern int CountIndepLogicalChannels(int *chans, int *nphys, int nlevels);

// Allocate and Free memory for the tree of lower-level channels which determine
// a logical channel.
extern int MemManageChannels(double *****channels, int *nphys, int *nencs, int nlevels, int decoder, int tofree);

// Allocate and free memory for the input channels structure in
// ComputeLogicalChannels(...).
extern void MemManageInputChannels(double ****inputchannels, int nphys, int nlogs, int decoder, int tofree);

// Allocate memory for the bins according to which logical (effective) channels
// must be averaged at intermediate decoding levels.
extern void AllocDecoderBins(struct simul_t *simul, int *nphys);

// Free memory allocated for decoding bins.
extern void FreeDecoderBins(struct simul_t *simul);

// Allocate memory for a 2D double array.
extern void AllocateDoubleArray2D(double **arr, int nrows, int ncols);

// Free memory for a 2D double array.
extern void FreeDoubleArray2D(double **arr, int nrows);

// Allocate memory for a 2D double complex array.
extern void AllocateDoubleComplexArray2D(double complex **arr, int nrows, int ncols);

// Free memory for a 2D double complex array.
extern void FreeDoubleComplexArray2D(double complex **arr, int nrows);

#endif /* MEMORY_H */
