#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include "constants.h"
#include "memory.h"
#include "printfuns.h"
#include "linalg.h"
#include "qecc.h"
#include "effective.h"
#include "benchmark.h"
// #include "logmetrics.h" // only for testing

void InitBenchOut(struct BenchOut *pbout, int nlevels, int nmetrics, int nlogs, int nbins, int nbreaks)
{
	// Initialize the memory allocated the to the benchmark output.
	pbout->logchans = malloc(sizeof(double) * (nlevels + 1) * nlogs * nlogs);
	pbout->chanvar = malloc(sizeof(double) * (nlevels + 1) * nlogs * nlogs);
	pbout->logerrs = malloc(sizeof(double) * (nlevels + 1) * nmetrics);
	pbout->logvars = malloc(sizeof(double) * (nlevels + 1) * nmetrics);
	pbout->bins = malloc(sizeof(int) * nmetrics * (nlevels + 1) * nbins * nbins);
	pbout->running = malloc(sizeof(double) * nmetrics * nbreaks);
}

void FreeBenchOut(struct BenchOut *pbout)
{
	// free the memory allocated to the Benchmark output.
	free(pbout->logchans);
	free(pbout->chanvar);
	free(pbout->logerrs);
	free(pbout->logvars);
	free(pbout->bins);
	free(pbout->running);
}

struct BenchOut Benchmark(int nlevels, int *nkd, int *SS, int *normalizer, double *normphases_real, double *normphases_imag, char *chname, int iscorr, double *physical, int rc, int nmetrics, char **metrics, int *decoders, int *dclookups, double *mpinfo, int *operators_LST, int hybrid, int *decoderbins, int *ndecoderbins, int frame, int nbreaks, long *stats, int nbins, int maxbin, int importance, double *refchan, double infidelity)
{
	/*
	Benchmark an error correcting scheme.
	Inputs:
		(a). Specifications of the error correcting code
			   1. number of concatenation layers: int levels
			   2. N,K,D for each code
					int **nkd : 2D array of integers
						where the i-th row gives the N,K,D values of the i-th code.
			   3. Stabilizer syndrome signs for each code
					 int ***SS: 3D array of integers
							where the i-th array gives the stabilizer syndrome signs for the i-th code.
			   4. Logical action each code
					 int ****normalizer: 4D array of integers
							where the i-th array gives the logical action for the i-th code.
			   5. Logical action phases for each code.
					 int ***normphases: 3D array where the i-th array gives the phases for the i-th code.
		(b). Error channel
			   1. Channel name: char *chname
			   2. Correlated noise flag: int iscorr
			   2. Physical noise map (as a process matrix)
					 double **physical: 2D double array
		(c). Metrics
			   1. Number of metrics to be computed: int nmetrics
			   2. Names of the metrics to be computed: char *metrics.
		(d). Specifications of the simulation
			   1. Type of syndrome sampling
					 int importance: string specifying the sampling algorithm.
			   2. Number of syndromes to be sampled at the top level, with breakpoints. long *stats: array of longs, where the i-th element gives the i-th breakpoint for the running average.
			   3. Number of syndrome metric bins: int nbins
			   4. Maximum order of magnitude for a bin: int maxbin
			   5. Quantum error correction frame: int frame
		(e). Decoding technique
			   1. Decoding technique, soft or hybrid: int hybrid: int hybrid
			   2. Channels at intermediate levels that need to be averaged: int **decoderbins
			   3. Number of distinct (bins) channels at each intermediate
	level: int *ndecoderbins
	*/
	
	/*
		iscorr
		physical
		decoders
		mpinfo
		importance
		infidelity
	*/
	struct constants_t *consts = malloc(sizeof(struct constants_t));
	InitConstants(consts);

	// Initialize the error correcting code structure.
	// printf("Quantum error correcting code with %d levels.\n", nlevels);
	struct qecc_t **qcode = malloc(sizeof(struct qecc_t *) * nlevels);
	int log, t, l, s, g, i, q, s_count = 0, ss_count = 0, normcount = 0, norm_phcount = 0, lst_count = 0, nparams = 0;
	for (l = 0; l < nlevels; l++)
	{
		// printf("l = %d\n", l);
		qcode[l] = malloc(sizeof(struct qecc_t));
		qcode[l]->N = nkd[3 * l];
		qcode[l]->K = nkd[3 * l + 1];
		qcode[l]->D = nkd[3 * l + 2];
		InitQECC(qcode[l]);
		for (s = 0; s < qcode[l]->nstabs; s++)
			for (g = 0; g < qcode[l]->nstabs; g++)
				(qcode[l]->projector)[s][g] = SS[ss_count + s * qcode[l]->nstabs + g];
		ss_count += qcode[l]->nstabs * qcode[l]->nstabs;

		// printf("ss_count = %d\n", ss_count);

		for (i = 0; i < qcode[l]->nlogs; i++)
			for (s = 0; s < qcode[l]->nstabs; s++)
				for (q = 0; q < qcode[l]->N; q++)
					(qcode[l]->action)[i][s][q] = normalizer[normcount + i * qcode[l]->nstabs * qcode[l]->N + s * qcode[l]->N + q];
		normcount += qcode[l]->nlogs * qcode[l]->nstabs * qcode[l]->N;

		// PrintIntArray2D((qcode[l]->action)[0], "S", qcode[l]->nstabs, qcode[l]->N);
		// PrintIntArray2D((qcode[l]->action)[1], "X S", qcode[l]->nstabs, qcode[l]->N);
		// PrintIntArray2D((qcode[l]->action)[2], "Y S", qcode[l]->nstabs, qcode[l]->N);
		// PrintIntArray2D((qcode[l]->action)[3], "Z S", qcode[l]->nstabs, qcode[l]->N);

		// printf("normcount = %d\n", normcount);

		for (i = 0; i < qcode[l]->nlogs; i++)
			for (s = 0; s < qcode[l]->nstabs; s++)
				(qcode[l]->phases)[i][s] = normphases_real[norm_phcount + i * qcode[l]->nstabs + s] + I * normphases_imag[norm_phcount + i * qcode[l]->nstabs + s];

		// printf("norm_phcount = %d\n", norm_phcount);

		if (decoders[l] == 1){
			// Lookup table for the minimum weight decoder
			for (s = 0; s < qcode[l]->nstabs; s++)
				(qcode[l]->dclookup)[s] = dclookups[s_count + s];
			s_count += qcode[l]->nstabs;
		}
		else;

		norm_phcount += qcode[l]->nlogs * qcode[l]->nstabs;

		for (log = 0; log < qcode[l]->nlogs; log ++)
			for (s = 0; s < qcode[l]->nstabs; s ++)
				for (t = 0; t < qcode[l]->nstabs; t ++)
					for (q = 0; q < qcode[l]->N; q ++)
						(qcode[l]->LST)[log][s][t][q] = operators_LST[lst_count + log * qcode[l]->nstabs * qcode[l]->nstabs * qcode[l]->N + s * qcode[l]->nstabs * qcode[l]->N + t * qcode[l]->N + q];
		lst_count += (int) pow(4, qcode[l]->N) * qcode[l]->N;
		
		// PrintIntArray2D((qcode[l]->TLS)[0][0], "T_0 L_0 S", qcode[l]->nstabs, qcode[l]->N);

		// printf("Code at level %d: N = %d, K = %d, D = %d.\n", l, qcode[l]->N, qcode[l]->K, qcode[l]->D);
	}

	// printf("QECC assigned.\n");

	// if (iscorr == 0)
	// 	PrintDoubleArray1D(physical, "physical channel", qcode[0]->nlogs * qcode[0]->nlogs);
	// else if (iscorr == 1)
	// 	PrintDoubleArray1D(physical, "physical channel", qcode[0]->nlogs * qcode[0]->nstabs);
	// else
	// 	PrintDoubleArray1D(physical, "physical channel", qcode[0]->N * qcode[0]->nlogs * qcode[0]->nlogs);

	// Record the number of physical qubits -- this is for setting the sizes of the decoder bins.
	int *nphys = malloc(sizeof(int) * nlevels);
	for (l = 0; l < nlevels; l++)
		nphys[l] = qcode[l]->N;
	int *chans = malloc(sizeof(int) * nlevels);
	CountIndepLogicalChannels(chans, nphys, nlevels);
	// Parameters that are specific to the Montecarlo simulations to estimate the logical error rate.

	// double *mpinfo_file, *physical_file;
	// double diff;

	struct simul_t **sims = malloc(sizeof(struct simul_t *) * (1 + (int)(decoders[0] == 2)));
	int m, j, c, chan_count = 0;
	for (s = 0; s <= (int) (decoders[0] == 2); s++)
	{
		sims[s] = malloc(sizeof(struct simul_t));
		sims[s]->nlevels = nlevels;
		sims[s]->nmetrics = nmetrics;
		sims[s]->importance = importance;
		sims[s]->hybrid = hybrid;
		sims[s]->nbins = nbins;
		sims[s]->maxbin = maxbin;
		sims[s]->nbreaks = nbreaks;
		sims[s]->nstats = stats[nbreaks - 1];
		if (s == 1)
			sims[s]->iscorr = 1;
		else
			sims[s]->iscorr = iscorr;

		// printf("Allocating simulation parameters for\ns = %d, nlevels = %d, nmetrics = %d, importance = %d, decoder = %d, nbins = %d, maxbin = %d, nstats = %ld, nbreaks = %d.\n", s, sims[s]->nlevels, sims[s]->nmetrics, sims[s]->importance, sims[s]->hybrid, sims[s]->nbins, sims[s]->maxbin, sims[s]->nstats, sims[s]->nbreaks);

		AllocSimParams(sims[s], qcode[0]->N, qcode[0]->K);

		// Type of decoding algorithm to be used.
		for (l = 0; l < nlevels; l++)
			(sims[s]->decoders)[l] = decoders[l];
		// Logical frame for Quantum error correction
		if (frame == 0)
			for (l = 0; l < nlevels; l++)
				(sims[s]->frames)[l] = 4;
		else if (frame == 1)
			for (l = 0; l < nlevels; l++)
				(sims[s]->frames)[l] = consts->nclifford;
		else
		{
			for (l = 0; l < nlevels - 1; l++)
				(sims[s]->frames)[l] = 4;
			(sims[s]->frames)[nlevels - 1] = consts->nclifford;
		}

		// printf("Allocating decoding bins for decoder %d.\n", sims[s]->hybrid);

		// Prescription for averaging channels at intermediate decoding levels
		if (sims[s]->hybrid > 0)
		{
			AllocDecoderBins(sims[s], nphys);
			chan_count = 0;
			for (l = 0; l < nlevels; l++)
			{
				for (c = 0; c < chans[l]; c++)
					(sims[s]->decbins)[l][c] = decoderbins[chan_count + c];
				chan_count += chans[l];
				(sims[s]->ndecbins)[l] = ndecoderbins[l];
			}
		}

		// Error model and metrics
		// printf("Loading error channel, iscorr = %d.\n", iscorr);
		sprintf(sims[s]->chname, "%s", chname);
		if (sims[s]->iscorr == 0)
			nparams = qcode[0]->nlogs * qcode[0]->nlogs;
		else if (sims[s]->iscorr == 1)
			nparams = qcode[0]->nlogs * qcode[0]->nstabs;
		else if (sims[s]->iscorr == 2)
			nparams = qcode[0]->N * qcode[0]->nlogs * qcode[0]->nlogs;
		else
			nparams = qcode[0]->nlogs * qcode[0]->nstabs * qcode[0]->nlogs * qcode[0]->nstabs;

		for (i = 0; i < nparams; i++)
		{
			if (s == 0)
				(sims[s]->physical)[i] = physical[i];
			else
				(sims[s]->physical)[i] = refchan[i];

			if (sims[s]->iscorr == 0)
				(sims[s]->logical)[0][i / (qcode[0]->nlogs)][i % (qcode[0]->nlogs)] = (sims[s]->physical)[i];
		}
		
		// printf("Loading %d metrics to be computed.\n", nmetrics);

		for (m = 0; m < nmetrics; m++)
		{
			sprintf((sims[s]->metricsToCompute)[m], "%s", metrics[m]);
			// (sims[s]->metricsToCompute)[m] = metrics[m];
			// printf("(sims[s]->metricsToCompute)[%d] = %s\n", m, (sims[s]->metricsToCompute)[m]);
		}

		// printf("Loading running statistics with %d breaks.\n", nbreaks);

		// Running average
		for (i = 0; i < nbreaks; i++)
			(sims[s]->runstats)[i] = stats[i];

		// Setting outlier syndrome probabilities.
		// Upper and lower limits for the probability of the outlier syndromes.
		// We will assume that the probability of the outlier syndromes is not more than p^d/2 and not less than 80% of p^d/2
		// For a physical noise process whose Pauli transfer matrix is G, we will define p = 0.5 + 0.5 * (4 - tr(G))/4.
		// Additionally, we want to make sure that 0.5 <= p <= 1. This is safe for the importance sampler since p ~ 0 will lead to an indefinite search in PowerSearch(...) in sampling.c.
		// We will follow the definition of infidelity in eq. 5.16 of https://arxiv.org/abs/1109.6887.pdf.
		// printf("infidelity = %g.\n", infidelity);
		if (infidelity == -1)
			infidelity = (4 - TraceFlattened(sims[s]->physical, qcode[0]->nlogs))/((double) 4);
		(sims[s]->outlierprobs)[1] = Max(0.4, infidelity);
		(sims[s]->outlierprobs)[0] = 0.80 * (sims[s]->outlierprobs)[1];

		// printf("infidelity = %g, Outlier probabilities lie in the range: [%g, %g].\n", infidelity, (sims[s]->outlierprobs)[0], (sims[s]->outlierprobs)[1]);

		// printf("Allocations complete for s = %d.\n", s);

		// Randomized compiling of quantum gates
		sims[s]->rc = rc;

		// Initial knowledge of Pauli error probabilities for a message passing decoder.
		// PrintIntArray1D((sims[s]->decoders), "Decoders", sims[s]->nlevels);
		if ((sims[s]->decoders)[0] == 3)
			for (i = 0; i < (int) pow(4, qcode[0]->N); i ++)
				(sims[s]->mpinfo)[i] = mpinfo[i];
		// mpinfo_file = malloc((int)pow(4, qcode[0]->N) * sizeof(double));
		// LoadDoubleArray1D(mpinfo_file, "./../input/debug_testing/mpinfo.txt", (int)pow(4, qcode[0]->N));
		// diff = 0;
		// for (i = 0; i < (int) pow(4, qcode[0]->N); i ++)
		// 	diff += fabs(mpinfo[i] - mpinfo_file[i]);
		// printf("Disparity in mpinfo = %g.\n", diff);
		// free(mpinfo_file);
		// physical_file = malloc(nparams * sizeof(double));
		// LoadDoubleArray1D(physical_file, "./../input/debug_testing/physical.txt", nparams);
		// diff = 0;
		// for (i = 0; i < nparams; i ++)
		// 	diff += fabs(physical[i] - physical_file[i]);
		// // printf("Disparity in physical = %g.\n", diff);
		// free(physical_file);
	}


	// printf("**************************************\n");
	// printf("INPUTS:\n");
	// printf("iscorr = %d\n", iscorr);
	// printf("nparams = %d\n", nparams);
	// PrintDoubleArray1D(physical, "Physical channel", nparams);
	// PrintIntArray1D(decoders, "Decoders", nlevels);
	// PrintDoubleArray1D(mpinfo, "Message passing information", (int)pow(4, qcode[0]->N));
	// printf("importance = %d\n", importance);
	// printf("infidelity = %.14f\n", infidelity);
	// printf("**************************************\n");

	// printf("Going to start Performance.\n");

	// ###################################

	Performance(qcode, sims, consts);

	// ###################################

	// printf("Loading outputs on to BenchOut.\n");

	int nlogs = qcode[0]->nlogs;
	struct BenchOut bout;
	InitBenchOut(&bout, sims[0]->nlevels, sims[0]->nmetrics, nlogs, sims[0]->nbins, (int)stats[0]);

	// Benchmark output
	// 1. Average effective channels for each concatenation level.
	// 2. Variance in the average effective channel for each concatenation level.
	// 3. For each metric
	// (a). Logical metric values for each concatenation level.
	// (b). Variance of logical error metric values for each concatenation level.
	// (c). Counts binned according to the syndrome probability and conditioned
	// logical error metric.
	// 4. Running average of average metric values for topmost level.

	// printf("Assigning output values.\n");

	for (l = 0; l < nlevels + 1; l++)
	{
		printf("l = %d\n", l);
		PrintDoubleArray2D((sims[0]->logical)[l], "logical channel", nlogs, nlogs);
		for (i = 0; i < nlogs; i++)
		{
			for (j = 0; j < nlogs; j++)
			{
				(bout.logchans)[l * (int)pow((double)nlogs, 2) + i * nlogs + j] = (sims[0]->logical)[l][i][j];
				(bout.chanvar)[l * (int)pow((double)nlogs, 2) + i * nlogs + j] = (sims[0]->variance)[l][sims[0]->nmetrics + i * nlogs + j];
			}
		}
	}
	// PrintDoubleArray1D((bout.logchans), "Logical channels", (nlevels + 1) * nlogs * nlogs);

	for (m = 0; m < nmetrics; m++)
	{
		for (l = 0; l < nlevels + 1; l++)
		{
			(bout.logerrs)[m * (nlevels + 1) + l] = (sims[0]->metricValues)[l][m];
			(bout.logvars)[m * (nlevels + 1) + l] = (sims[0]->variance)[l][m];
			for (i = 0; i < nbins; i++)
				for (j = 0; j < nbins; j++)
					(bout.bins)[m * (nlevels + 1) * (int)pow((double)nbins, 2) + l * (int)pow((double)nbins, 2) + i * nbins + j] = (sims[0]->bins)[l][m][i][j];
		}
		for (i = 0; i < nbreaks; i++)
			(bout.running)[m * nbreaks + i] = (sims[0]->runavg)[m][i + 1];
	}
	PrintDoubleArray1D((bout.logerrs), "Metric values", nmetrics * (nlevels + 1));
	// PrintDoubleArray1D((bout.running), "Running average", nmetrics * nbreaks);

	// ######################################
	// printf("Freeing memory.\n");

	// Free memory
	free(nphys);
	free(chans);
	// printf("Freeing %d simulation structures.\n", 1 + (int)(importance == 2));
	for (s = 0; s < 1 + (int)(decoders[0] == 2); s++)
	{
		FreeSimParams(sims[s], qcode[0]->N, qcode[0]->K);
		if (sims[s]->hybrid > 0)
			FreeDecoderBins(sims[s]);
		free(sims[s]);
	}
	free(sims);
	// printf("Freeing %d qcode structures.\n", nlevels);
	for (l = 0; l < nlevels; l++)
	{
		FreeQECC(qcode[l]);
		free(qcode[l]);
	}
	free(qcode);
	FreeConstants(consts);

	// printf("Benchmark complete!\n");
	return bout;
}
