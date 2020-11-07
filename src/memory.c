#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "constants.h"
#include "printfuns.h"
#include "rand.h"
#include "linalg.h"
#include "memory.h"

void AllocateDoubleArray2D(double **arr, int nrows, int ncols){
	// Allocate memory for a 2D double array.
	int r, c;
	arr = malloc(sizeof(double *) * nrows);
	for (r = 0; r < nrows; r ++){
		arr[r] = malloc(sizeof(double) * ncols);
		for (c = 0; c < ncols; c ++){
			arr[r][c] = 0;
		}
	}
}

void FreeDoubleArray2D(double **arr, int nrows){
	// Free memory for a 2D double array.
	int r;
	for (r = 0; r < nrows; r ++){
		free(arr[r]);
	}
	free(arr);
}

void AllocateDoubleComplexArray2D(double complex **arr, int nrows, int ncols){
	// Allocate memory for a 2D double complex array.
	int r, c;
	arr = malloc(sizeof(double complex *) * nrows);
	for (r = 0; r < nrows; r ++){
		arr[r] = malloc(sizeof(double complex) * ncols);
		for (c = 0; c < ncols; c ++){
			arr[r][c] = 0 + 0 * I;
		}
	}
}

void FreeDoubleComplexArray2D(double complex **arr, int nrows){
	// Free memory for a 2D double complex array.
	int r;
	for (r = 0; r < nrows; r ++){
		free(arr[r]);
	}
	free(arr);
}


void AllocSimParamsQECC(struct simul_t *simul, int nphys, int nenc)
{
	// Allocate memory for parameters of the simulation structure that depend upon
	// the QECC used. This memeory will be reallocated everytime there is a new
	// QECC, i.e, at a new concatenation level. These parameters are virtchan,
	// logical, syndprobs, cumulative, levelZeroSynds, levelZeroCumulative,
	// levelZeroImpDist, levelZeroImpCumul, process, corrections, effprocess,
	// effective, levelZeroChannels
	int q, i, nstabs = (int)pow(2, ((double)(nphys - nenc))), nlogs = (int)pow(4, (double)nenc);
	// printf("Function: AllocSimParamsQECC, nstabs = %d, nlogs = %d.\n", nstabs, nlogs);
	simul->virtchan = malloc(nphys * sizeof(double **));
	simul->pauli_probs = malloc(nphys * sizeof(double *));
	for (q = 0; q < nphys; q++)
	{
		(simul->virtchan)[q] = malloc(nlogs * sizeof(double *));
		(simul->pauli_probs)[q] = malloc(nlogs * sizeof(double));
		for (i = 0; i < nlogs; i++)
			(simul->virtchan)[q][i] = malloc(nlogs * sizeof(double));
	}
	// printf("_/ virtchan\n");
	// simul->rcpauli = malloc((nstabs) * sizeof(int));
	// int pure = RandomRangeInt(0, nstabs);
	// // printf("Applying RC by compiling with the pure error %d.\n", pure);
	// for (i = 0; i < nstabs; i ++)
	//  	(simul->rcpauli)[i] = BinaryDot(i, pure);
	// Syndrome sampling
	(simul->syndprobs) = malloc(nstabs * sizeof(double));
	(simul->cumulative) = malloc(nstabs * sizeof(double));
	int s;
	for (s = 0; s < nstabs; s++)
	{
		(simul->syndprobs)[s] = 0.0;
		(simul->cumulative)[s] = 0.0;
	}
	// printf("_/ syndprobs, cumulative\n");
	// Quantum error correction
	simul->process = malloc(nlogs * sizeof(double ***));
	int j;
	for (i = 0; i < nlogs; i++)
	{
		(simul->process)[i] = malloc(nlogs * sizeof(double **));
		for (j = 0; j < nlogs; j++)
		{
			(simul->process)[i][j] = malloc(nstabs * sizeof(double *));
			for (s = 0; s < nstabs; s++)
				(simul->process)[i][j][s] = malloc(nstabs * sizeof(double));
		}
	}
	// printf("_/ process\n");
	(simul->corrections) = malloc(nstabs * sizeof(int));
	// printf("_/ corrections\n");
	(simul->effprocess) = malloc(nstabs * sizeof(double **));
	(simul->cosetprobs) = malloc(nstabs * sizeof(double *));
	for (s = 0; s < nstabs; s++)
	{
		(simul->effprocess)[s] = malloc(nlogs * sizeof(double *));
		(simul->cosetprobs)[s] = malloc(nlogs * sizeof(double));
		for (i = 0; i < nlogs; i++)
		{
			(simul->cosetprobs)[s][i] = 0;
			(simul->effprocess)[s][i] = malloc(nlogs * sizeof(double));
		}
	}
	// printf("_/ effective, effprocess\n");
}

void AllocDecoderBins(struct simul_t *simul, int *nphys)
{
	// Channels to be averaged over, at intermediate levels of the decoding
	// process. Channels will be averaged in the Coarsegrain(...) method of
	// effective.c. The channels to be averaged all have the same index in this
	// array.
	int *chans = malloc(sizeof(int) * simul->nlevels);
	CountIndepLogicalChannels(chans, nphys, simul->nlevels);
	simul->decbins = malloc(sizeof(int *) * simul->nlevels);
	simul->ndecbins = malloc(sizeof(int) * simul->nlevels);
	int l;
	for (l = 0; l < simul->nlevels; l++)
		(simul->decbins)[l] = malloc(sizeof(int) * chans[l]);
	free(chans);
}

void FreeDecoderBins(struct simul_t *simul)
{
	// Free the memory allocated to decoder bins
	int l;
	for (l = 0; l < simul->nlevels; l++)
		free((simul->decbins)[l]);
	free(simul->decbins);
	free(simul->ndecbins);
}

void AllocSimParams(struct simul_t *simul, int nphys, int nenc)
{
	// Initialize the elements that pertain to the montecarlo simulation of
	// channels. Physical channels.
	int i, nstabs = (int)pow(2, ((double)(nphys - nenc))), nlogs = (int)pow(4, (double)nenc);
	// printf("Function: AllocSimParams, nstabs = %d, nlogs = %d\n", nstabs, nlogs);
	simul->chname = malloc(100 * sizeof(char));
	int nparams;
	if (simul->iscorr == 0)
		nparams = nlogs * nlogs;
	else if (simul->iscorr == 1)
		nparams = nlogs * nstabs;
	else if (simul->iscorr == 2)
		nparams = nphys * nlogs * nlogs;
	else
		nparams = nlogs * nstabs * nlogs * nstabs;
	simul->physical = malloc(nparams * sizeof(double));
	for (i = 0; i < nparams; i++)
		simul->physical[i] = 0;

	// Initial guess for the message passing decoder.
	simul->mpinfo = malloc(sizeof(double) * (int) pow(4, nphys));

	// Metrics to be computed at every level.
	int m, l;
	simul->metricsToCompute = malloc(simul->nmetrics * sizeof(char *));
	for (m = 0; m < simul->nmetrics; m++)
	{
		(simul->metricsToCompute)[m] = malloc(100 * sizeof(char));
		// printf("metricsToCompute[%d] = %s.\n", m, (simul->metricsToCompute)[m]);
	}
	simul->metricValues = malloc((simul->nlevels + 1) * sizeof(double *));
	for (l = 0; l < simul->nlevels + 1; l++)
	{
		(simul->metricValues)[l] = malloc(simul->nmetrics * sizeof(double));
		for (m = 0; m < simul->nmetrics; m++)
			(simul->metricValues)[l][m] = 0;
	}

	// Average logical channel at top level.
	simul->logical = malloc((simul->nlevels + 1) * sizeof(double **));
	int j;
	for (l = 0; l < simul->nlevels + 1; l++)
	{
		(simul->logical)[l] = malloc(nlogs * sizeof(double *));
		for (i = 0; i < nlogs; i++)
		{
			(simul->logical)[l][i] = malloc(nlogs * sizeof(double));
			for (j = 0; j < nlogs; j++)
				(simul->logical)[l][i][j] = 0;
		}
	}

	// Syndrome sampling.
	simul->outlierprobs = malloc(2 * sizeof(double));

	simul->levelOneSynds = malloc(nstabs * sizeof(double));
	simul->levelOneImpDist = malloc(nstabs * sizeof(double));
	simul->levelOneCumul = malloc(nstabs * sizeof(double));
	simul->levelOneImpCumul = malloc(nstabs * sizeof(double));
	simul->statsperlevel = malloc((simul->nlevels + 1) * sizeof(long));
	for (l = 0; l < simul->nlevels + 1; l++)
		(simul->statsperlevel)[l] = 0;

	// Syndrome-metric bins.
	simul->bins = malloc((simul->nlevels + 1) * sizeof(int ***));
	for (l = 0; l < simul->nlevels + 1; l++)
	{
		(simul->bins)[l] = malloc((simul->nmetrics) * sizeof(int **));
		for (m = 0; m < simul->nmetrics; m++)
		{
			(simul->bins)[l][m] = malloc(simul->nbins * sizeof(int *));
			for (i = 0; i < simul->nbins; i++)
			{
				(simul->bins)[l][m][i] = malloc(simul->nbins * sizeof(int));
				for (j = 0; j < simul->nbins; j++)
					(simul->bins)[l][m][i][j] = 0;
			}
		}
	}

	// Variance measures.
	simul->sumsq = malloc((simul->nlevels + 1) * sizeof(double *));
	simul->variance = malloc((simul->nlevels + 1) * sizeof(double *));
	for (l = 0; l < simul->nlevels + 1; l++)
	{
		(simul->sumsq)[l] = malloc((simul->nmetrics + nlogs * nlogs) * sizeof(double));
		(simul->variance)[l] = malloc((simul->nmetrics + nlogs * nlogs) * sizeof(double));
		for (i = 0; i < simul->nmetrics + nlogs * nlogs; i++)
		{
			(simul->sumsq)[l][i] = 0;
			(simul->variance)[l][i] = 0;
		}
	}

	// fprintf(stderr, "_/ sumsq, variance.\n");
	// PrintDoubleArray2D(simul[0].sumsq, "sumsq", simul[0].nlevels + 1,
	// simul[0].nmetrics + nlogs * nlogs).

	// Running statistics
	simul->runstats = malloc(sizeof(long) * simul->nbreaks);
	simul->runavg = malloc(sizeof(double *) * simul->nmetrics);
	for (m = 0; m < simul->nmetrics; m++)
	{
		(simul->runavg)[m] = malloc(sizeof(double) * (simul->nbreaks + 1));
		for (i = 0; i < simul->nbreaks + 1; i++)
			(simul->runavg)[m][i] = 0;
	}

	// Quantum error correction.
	simul->levelOneChannels = malloc(nstabs * sizeof(double **));
	simul->levelOneCosets = malloc(nstabs * sizeof(double *));
	int s;
	for (s = 0; s < nstabs; s++)
	{
		(simul->levelOneChannels)[s] = malloc(nlogs * sizeof(double *));
		(simul->levelOneCosets)[s] = malloc(nlogs * sizeof(double));
		for (i = 0; i < nlogs; i++)
			(simul->levelOneChannels)[s][i] = malloc(nlogs * sizeof(double));
	}
	simul->frames = malloc(simul->nlevels * sizeof(int));
	// printf("_/ levelOneChannels, frames.\n");

	// Specify the type of decoder being used: 0 for maximum likelihood and 1 for minimum weight.
	simul->decoders = malloc(simul->nlevels * sizeof(int));
	for (i = 0; i < simul->nlevels; i++)
		(simul->decoders)[i] = 0;
}

void FreeSimParamsQECC(struct simul_t *simul, int nphys, int nenc)
{
	// Free the memory allocated to simulation parameters that depend on the QECC.
	// These parameters are
	// virtchan, logical, syndprobs, cumulative, levelZeroSynds,
	// levelZeroCumulative, levelZeroImpDist, levelZeroImpCumul, process,
	// corrections, effprocess, effective, levelZeroChannels.
	// printf("FreeSimParamsQECC: nphys = %d, nstabs = %d, nlogs = %d\n", nphys,
	// nstabs, nlogs); Logical channels at intermediate levels.
	int q, i, nstabs = (int)pow(2, ((double)(nphys - nenc))), nlogs = (int)pow(4, (double)nenc);
	for (q = 0; q < nphys; q++)
	{
		for (i = 0; i < nlogs; i++)
			free((simul->virtchan)[q][i]);
		free((simul->virtchan)[q]);
		free((simul->pauli_probs)[q]);
	}
	free(simul->virtchan);
	free(simul->pauli_probs);
	// free(simul->rcpauli);
	// printf("_/ virtchan\n");
	// Quantum error correction.
	int j, s;
	for (i = 0; i < nlogs; i++)
	{
		for (j = 0; j < nlogs; j++)
		{
			for (s = 0; s < nstabs; s++)
				free((simul->process)[i][j][s]);
			free((simul->process)[i][j]);
		}
		free((simul->process)[i]);
	}
	free(simul->process);
	free(simul->corrections);
	// printf("_/ process, corrections\n");
	for (s = 0; s < nstabs; s++)
	{
		for (i = 0; i < nlogs; i++)
		{
			free((simul->effprocess)[s][i]);
		}
		free(simul->effprocess[s]);
		free(simul->cosetprobs[s]);
	}
	free(simul->effprocess);
	free(simul->cosetprobs);
	// printf("_/ effective, effprocess\n");
	// Syndrome sampling.
	free(simul->syndprobs);
	free(simul->cumulative);
	// printf("_/ syndprobs, cumulative\n");
	// printf("Done\n")
}

void FreeSimParams(struct simul_t *simul, int nphys, int nenc)
{
	// Free memory allocated to the simulation structure.
	// Physical channels.
	int i, nstabs = (int)pow(2, ((double)(nphys - nenc))), nlogs = (int)pow(4, (double)nenc);
	// printf("Function: FreeSimParams with nstabs = %d, nlogs = %d\n", nstabs, nlogs);
	free(simul->physical);
	free(simul->mpinfo);
	free(simul->chname);
	// printf("Freed physical and chname.\n");
	// Metrics to be computed at logical levels.
	int m;
	for (m = 0; m < simul->nmetrics; m++)
		free((simul->metricsToCompute)[m]);
	free(simul->metricsToCompute);
	int l;
	for (l = 0; l < 1 + simul->nlevels; l++)
		free((simul->metricValues)[l]);
	free(simul->metricValues);
	// printf("Freed metricsToCompute and metricValues.\n");
	// Average logical channel at the top level.
	for (l = 0; l < 1 + simul->nlevels; l++)
	{
		for (i = 0; i < nlogs; i++)
			free((simul->logical)[l][i]);
		free((simul->logical)[l]);
	}
	free(simul->logical);
	// printf("Freed logical.\n");
	// Syndrome sampling.
	free(simul->outlierprobs);
	free(simul->statsperlevel);
	free(simul->levelOneSynds);
	free(simul->levelOneCumul);
	free(simul->levelOneImpDist);
	free(simul->levelOneImpCumul);
	// printf("Freed statsperlevel, levelOneSynds, levelOneCumul, levelOneImpDist,
	// levelOneCumul.\n"); Variance measure.
	for (l = 0; l < 1 + simul->nlevels; l++)
		free((simul->sumsq)[l]);
	free(simul->sumsq);
	for (l = 0; l < 1 + simul->nlevels; l++)
		free((simul->variance)[l]);
	free(simul->variance);
	// printf("Freed sumsq, variance.\n");
	// Running statistics
	free(simul->runstats);
	for (m = 0; m < simul->nmetrics; m++)
		free((simul->runavg)[m]);
	free(simul->runavg);
	// printf("Freed runstats and runavg.\n");
	// Syndrome metric bins.
	for (l = 0; l < 1 + simul->nlevels; l++)
	{
		for (m = 0; m < simul->nmetrics; m++)
		{
			for (i = 0; i < simul->nbins; i++)
				free((simul->bins)[l][m][i]);
			free((simul->bins)[l][m]);
		}
		free((simul->bins)[l]);
	}
	free(simul->bins);
	// printf("Freed bins.\n");
	// Quantum error correction.
	int s;
	for (s = 0; s < nstabs; s++)
	{
		// printf("s = %d:\n", s);
		// PrintDoubleArray2D((simul->levelOneChannels)[s], "Level 1 channel", nlogs, nlogs);
		for (i = 0; i < nlogs; i++)
			free((simul->levelOneChannels)[s][i]);
		free((simul->levelOneChannels)[s]);
		free((simul->levelOneCosets)[s]);
	}
	free(simul->levelOneChannels);
	free(simul->levelOneCosets);
	// printf("Freed levelOneChannels.\n");
	free(simul->frames);
	free(simul->decoders);
	// printf("Freed frames.\n");
}

int CountIndepLogicalChannels(int *chans, int *nphys, int nlevels)
{
	// Determine the number of independent logical channels at every level, that
	// determine the logical channels of higher levels.
	// printf("Function: CountIndepLogicalChannels\n")
	chans[nlevels-1] = 1;
	int l;
	for (l = nlevels - 2; l >= 0; l--)
		chans[l] = nphys[l] * chans[l + 1];
	return chans[0];
}

int MemManageChannels(double *****channels, int *nphys, int *nencs, int nlevels, int decoder, int tofree)
{
	// Allocate and Free memory for the tree of lower-level channels which
	// determine a logical channel. printf("Function: MemManageChannels, tofree =
	// %d\n", tofree);
	int c, j, s, l, nchans, nlogs;
	int *chans = malloc(sizeof(int) * (nlevels));
	nchans = CountIndepLogicalChannels(chans, nphys, nlevels);
	// PrintIntArray1D(chans, "Channels", nlevels);
	// printf("Creating %d channels.\n", nchans);
	if (tofree == 0)
	{
		// Allocate memory.
		for (l = 0; l < nlevels; l++)
		{
			channels[l] = malloc(sizeof(double ***) * chans[l]);
			nlogs = (int)pow(4, (double)nencs[l]);
			for (c = 0; c < chans[l]; c++)
			{
				channels[l][c] = malloc(sizeof(double **) * (1 + (int)(decoder == 2)));
				for (s = 0; s < 1 + (int)(decoder == 2); s++)
				{
					channels[l][c][s] = malloc(sizeof(double *) * (2 + nlogs));
					for (j = 0; j < nlogs; j++)
						channels[l][c][s][j] = malloc(sizeof(double) * nlogs);
					channels[l][c][s][nlogs] = malloc(sizeof(double) * 3);
					channels[l][c][s][nlogs + 1] = malloc(sizeof(double) * nlogs);
				}
			}
		}
	}
	else
	{
		// Free memory
		for (l = 0; l < nlevels; l++)
		{
			nlogs = (int)pow(4, (double)nencs[l]);
			for (c = 0; c < chans[l]; c++)
			{
				for (s = 0; s < 1 + (int)(decoder == 2); s++)
				{
					for (j = 0; j < 2 + nlogs; j++)
						free(channels[l][c][s][j]);
					free(channels[l][c][s]);
				}
				free(channels[l][c]);
			}
			free(channels[l]);
		}
	}
	free(chans);
	return nchans;
}

void MemManageInputChannels(double ****inputchannels, int nphys, int nlogs, int decoder, int tofree)
{
	// Allocate and free memory for the input channels structure in
	// ComputeLogicalChannels(...). printf("Function: MemManageInputChannels,
	// tofree = %d, nlogs = %d.\n", tofree, nlogs);
	int i, q, s;
	if (tofree == 0)
	{
		// Initialize the space required for the input channels.
		for (q = 0; q < nphys; q++)
		{
			inputchannels[q] = malloc(sizeof(double **) * (1 + (int)(decoder == 2)));
			for (s = 0; s < 1 + (int)(decoder == 2); s++)
			{
				inputchannels[q][s] = malloc(sizeof(double *) * (nlogs + 1));
				for (i = 0; i < nlogs; i++)
					inputchannels[q][s][i] = malloc(sizeof(double) * nlogs);
				inputchannels[q][s][nlogs] = malloc(sizeof(double) * 2);
			}
		}
	}
	else
	{
		// Free memory
		for (q = 0; q < nphys; q++)
		{
			for (s = 0; s < 1 + (int)(decoder == 2); s++)
			{
				for (i = 0; i < 1 + nlogs; i++)
					free(inputchannels[q][s][i]);
				free(inputchannels[q][s]);
			}
			free(inputchannels[q]);
		}
	}
	// printf("Done MemManageInputChannels.\n");
}
