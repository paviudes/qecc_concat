#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "mt19937/mt19937ar.h"
#include "constants.h"
#include "printfuns.h"
#include "memory.h"
#include "qecc.h"
#include "checks.h"
#include "logmetrics.h"
#include "sampling.h"
#include "hybrid.h"
#include "decode.h"
#include "effective.h"

int IsElement(long *arr, int size, long item) {
	// Determine if an item is present in an array
	// printf("Testing if %d is present in the following array of size %d.\n",
	// item, arr[0]); PrintIntArray1D(arr, "array", arr[0] + 1);
	int i;
	for (i = 0; i < size; i++)
		if (arr[i] == item)
			return 1;
	return 0;
}

int GetBinPosition(double number, int nbins, int maxbin) {
	// Place a number into one of the bins depending on the order of magnitude of
	// the number. If a number is of the order of magnitude of 10^-i, then return
	// i. To find the order of magnitude we will take the negative log and bin it
	// from 0 to 20.
	int binindex = (int)((-1) * log10(number) / ((double)maxbin) * ((double)nbins));
	if (binindex < 0)
		binindex = 0;
	if (binindex > (nbins - 1))
		binindex = (nbins - 1);
	return binindex;
}

int PickFromInterval(int low, int high, int *arr, int *found, int alloc) {
	// Find the elements in the given (sorted) array that are in a specified
	// range. The found elements must be recorded in the "found" array. The first
	// elements of the sorted and found arrays indicate their respective sizes.
	// printf("Find elements in the array\n");
	// PrintIntArray1D(arr, "array", arr[0] + 1);
	// printf("that are between %d and %d.\n", low, high);
	int i = 0, nelems = 0;
	for (i = 0; i < arr[0]; i++) {
		if ((arr[i + 1] >= low) && (arr[i + 1] <= high)) {
			nelems = nelems + 1;
			if (alloc == 1)
				found[nelems] = arr[i + 1];
		}
	}
	if (alloc == 1)
		found[0] = nelems;
	return nelems;
}

void UpdateMetrics(int level, double bias, double history, int isfinal, struct qecc_t *qcode, struct simul_t *sim, struct constants_t *consts) {
	// Compute metrics for all the effective channels and update the average value
	// of the metrics. Metric values.
	// printf("Updating metrics.\n");
	double *metvals = malloc(sim->nmetrics * sizeof(double));
	double *avg = malloc((sim->nmetrics + qcode->nlogs * qcode->nlogs) * sizeof(double));
	int m;
	for (m = 0; m < sim->nmetrics + qcode->nlogs * qcode->nlogs; m++)
		avg[m] = 0;

	int s, i, j;
	if (isfinal == 0) {
		for (s = 0; s < qcode->nstabs; s++) {
			if ((sim->syndprobs)[s] > consts->atol) {
				// Compute metrics.
				ComputeMetrics(metvals, sim->nmetrics, sim->metricsToCompute, sim->effprocess[s], sim->chname, consts);
				// if (s < 2) {
				// 	printf("s = %d\n", s);
				// 	PrintDoubleArray1D(metvals, "metric values for s", sim->nmetrics);
				// }
				for (m = 0; m < sim->nmetrics; m++)
					avg[m] += metvals[m] * (sim->syndprobs)[s];
				// Compute average channel.
				for (i = 0; i < qcode->nlogs; i++)
					for (j = 0; j < qcode->nlogs; j++)
						avg[sim->nmetrics + i * qcode->nlogs + j] += (sim->effprocess)[s][i][j] * (sim->syndprobs)[s];
			}
		}

		// PrintDoubleArray1D(avg, "Avg Metric values and channels", sim->nmetrics + qcode->nlogs * qcode->nlogs);
		// Average of metrics.
		for (m = 0; m < sim->nmetrics; m++) {
			(sim->metricValues)[level + 1][m] += bias * avg[m];
			(sim->sumsq)[level + 1][m] += pow(bias * avg[m], 2);
		}
		// PrintDoubleArray1D((sim->metricValues)[level + 1], "(sim->metricValues)[level + 1]", sim->nmetrics);
		// printf("Bias = %g.\n", bias);
		// printf("Populating average logical channel.\n");
		for (i = 0; i < qcode->nlogs; i++) {
			for (j = 0; j < qcode->nlogs; j++) {
				(sim->logical)[level + 1][i][j] += bias * avg[sim->nmetrics + i * qcode->nlogs + j];
				(sim->sumsq)[level + 1][sim->nmetrics + i * qcode->nlogs + j] += pow(bias * avg[sim->nmetrics + i * qcode->nlogs + j], 2);
			}
		}
		// printf("The trace-preserving part of the logical channel at level %d is: %g and the bias is: %g.\n", level + 1, avg[sim->nmetrics], bias);

		// printf("Syndrome metric binning.\n");

		// Syndrome-metric binning.
		int mbin, sbin = GetBinPosition(fabsl(history), sim->nbins, sim->maxbin);
		for (m = 0; m < sim->nmetrics; m++) {
			mbin = GetBinPosition(fabsl(avg[m]), sim->nbins, sim->maxbin);
			(sim->bins)[level + 1][m][sbin][mbin] ++;
		}
		// Update the number of statistics done for the level.
		(sim->statsperlevel)[level + 1] ++;
	}
	else{
		printf("(sim->statsperlevel)[%d] = %ld.\n", level + 1, (sim->statsperlevel)[level + 1]);
		PrintDoubleArray1D(sim->metricValues[level + 1], "Metrics", sim->nmetrics);
		// After all the simulations are done, the average metrics are to be computed by diving the metricValues by the total number of statistics done for that level.
		for (m = 0; m < sim->nmetrics; m++) {
			(sim->metricValues)[level + 1][m] /= ((double)((sim->statsperlevel)[level + 1]));
			(sim->variance)[level + 1][m] = 1 / ((double)((sim->statsperlevel)[level + 1] * ((sim->statsperlevel)[level + 1] - 1))) * sim->sumsq[level + 1][m] - pow((sim->metricValues)[level + 1][m], 2);
		}
		for (i = 0; i < qcode->nlogs; i++) {
			for (j = 0; j < qcode->nlogs; j++) {
				(sim->logical)[level + 1][i][j] /= ((double)((sim->statsperlevel)[level + 1]));
				(sim->variance)[level + 1][sim->nmetrics + i * qcode->nlogs + j] = 1 / ((double)((sim->statsperlevel)[level + 1] * ((sim->statsperlevel)[level + 1] - 1))) * sim->sumsq[level + 1][sim->nmetrics + i * qcode->nlogs + j] - pow((sim->logical)[level + 1][i][j], 2);
			}
		}
		// printf("The trace-preserving part of the average logical channel at level %d is: %g.\n", level + 1, (sim->logical)[level + 1][0][0]);
	}
	// Free memory.
	free(metvals);
	free(avg);
	// printf("Updated metrics.\n");
}

void ComputeLevelZeroMetrics(struct simul_t *sim, int nqubits, int nlogs, struct constants_t *consts){
	// Compute the level-0 (physical) metrics.
	double **ptm, *chanmets;
	int i, j, q;

	// printf("nlogs = %d, nqubits = %d\n", nlogs, nqubits);

	if (sim->iscorr == 0){
		// ProcessToChoi((sim->logical)[0], nlogs, choi, consts->pauli);
		ComputeMetrics((sim->metricValues)[0], sim->nmetrics, sim->metricsToCompute, (sim->logical)[0], sim->chname, consts);
	}
	else if (sim->iscorr == 2){
		// Compute the level 0 metrics for each qubit's physical channel and average them.
		// We are assuming that the error metrics under with tensor products.
		// Extract the single qubit maps from the physical channels.
		// printf("Computing Level 0 metrics\n");
		// Fidelity is multiplicative whereas the other error-metrics are additive.
		// Infidelity is 1 - the product of fidelities
		for (i = 0; i < sim->nmetrics; i ++)
			if (strncmp((sim->metricsToCompute)[i], "infid", 5) == 0)
				(sim->metricValues)[0][i] = 1;

		// PrintDoubleArray1D((sim->metricValues)[0], "Level-0 metrics", sim->nmetrics);
		ptm = malloc(sizeof(double *) * nlogs);
		for (i = 0; i < nlogs; i++)
			ptm[i] = malloc(sizeof(double) * nlogs);

		chanmets = malloc(sizeof(double) * nqubits);
		for (q = 0; q < nqubits; q ++){
			// Extract the physical channel for the qubit q.
			for (i = 0; i < nlogs; i ++){
				for (j = 0; j < nlogs; j ++){
					ptm[i][j] = (sim->physical)[q * nlogs * nlogs + i * nlogs + j];
				}
			}
			// PrintDoubleArray2D(ptm, "Pauli transfer matrix", nlogs, nlogs);

			// Compute metrics for the q-th qubit's channel.
			for (i = 0; i < sim->nmetrics; i ++)
				chanmets[i] = 0;
			ComputeMetrics(chanmets, sim->nmetrics, sim->metricsToCompute, ptm, sim->chname, consts);
			// printf("qubit: %d\n", q);
			// PrintDoubleArray1D(chanmets, "Physical metrics", sim->nmetrics);
			// Add the metrics across qubits for all except fidelity.
			// For fidelity, we need to multiply the values.
			for (i = 0; i < sim->nmetrics; i ++){
				// printf("Metric: %s\n", sim->metricsToCompute[i]);
				if (strncmp(sim->metricsToCompute[i], "infid", 5) == 0)
					(sim->metricValues)[0][i] *= (1 - chanmets[i]);
				else
					(sim->metricValues)[0][i] += chanmets[i];
			}
		}
		// Infidelity is 1 - the product of fidelities
		for (i = 0; i < sim->nmetrics; i ++)
			if (strncmp(sim->metricsToCompute[i], "infid", 5) == 0)
				(sim->metricValues)[0][i] = 1 - (sim->metricValues)[0][i];
		// Free memory
		for (i = 0; i < nlogs; i++)
			free(ptm[i]);
		free(ptm);
		free(chanmets);
	}
	else;
	// PrintDoubleArray1D((sim->metricValues)[0], "Level-0 metrics", sim->nmetrics);
}

void ComputeLevelOneChannels(struct simul_t *sim, struct qecc_t *qcode, struct constants_t *consts, int decoder, int copy_lookup) {
	// Compute the effective channels and syndrome probabilities for all level-1
	// syndromes. The pre-computation is to avoid re-computing level-1 syndromes
	// for every new top-level syndrome. Load the physical channels on to the
	// simulation structure and perform qcode
	// printf("Function: ComputeLevelOneChannels\n");
	int q, i, j, isPauli = 1;
	// PrintIntArray1D((sim->decoders), "Decoders", sim->nlevels);
	// printf("Allocating resources level one for %d logical and %d physical qubits\n", qcode->K, qcode->N);
	AllocSimParamsQECC(sim, qcode->N, qcode->K);

	if ((sim->decoders)[0] == 3)
		ComputeCosetProbsLevelOne(sim->mpinfo, qcode->nlogs, qcode->nstabs, sim->cosetprobs);

	if ((sim->iscorr == 0) || (sim->iscorr == 2)) {
		for (q = 0; q < qcode->N; q++) {
			for (i = 0; i < qcode->nlogs; i++){
				for (j = 0; j < qcode->nlogs; j++){
					if (sim->iscorr == 0)
						(sim->virtchan)[q][i][j] = (sim->physical)[i * qcode->nlogs + j];
					else
						(sim->virtchan)[q][i][j] = (sim->physical)[q * qcode->nlogs * qcode->nlogs + i * qcode->nlogs + j];
				}
			}
			if (isPauli > 0)
				isPauli = isPauli * IsDiagonal((sim->virtchan)[q], qcode->nlogs);
			// printf("Qubit %d.\n", q);
			// PrintDoubleArrayDiag((sim->virtchan)[q], "Level 0 channel", qcode->nlogs);
		}

		// printf("Loaded virtual channels, isPauli = %d.\n", isPauli);
		SingleShotErrorCorrection(isPauli, sim->iscorr, decoder, (sim->frames)[0], qcode, sim, consts, 1);
	}
	else{
		if (sim->iscorr == 1)
			isPauli = 1;
		else
			isPauli = 0;
		// This runs for iscorr = 1 as well as iscorr = 3.
		// Correlated channel -- the physical channel contains the full process matrix
		// printf("Simulating a correlated channel: iscorr = %d, using decoder %d.\n", sim->iscorr, decoder);
		SetFullProcessMatrix(qcode, sim, sim->physical, isPauli);
		// printf("Running SingleShotErrorCorrection.\n");
		SingleShotErrorCorrection(isPauli, sim->iscorr, decoder, (sim->frames)[0], qcode, sim, consts, 1);
	}

	// printf("Completed SingleShotErrorCorrection.\n");

	UpdateMetrics(0, 1, 1, 0, qcode, sim, consts);
	int s;

	for (s = 0; s < qcode->nstabs; s++) {
		(sim->levelOneSynds)[s] = (sim->syndprobs)[s];
		for (i = 0; i < qcode->nlogs; i++){
			for (j = 0; j < qcode->nlogs; j++)
				(sim->levelOneChannels)[s][i][j] = (sim->effprocess)[s][i][j];
			(sim->levelOneCosets)[s][i] = (sim->cosetprobs)[s][i];
		}
		// printf("s = %d: P(s) = %g,\n", s, (sim->levelOneSynds)[s]);
		// PrintDoubleArray1D((sim->levelOneCosets)[s], "Level 1 coset probability", qcode->nlogs);
		// PrintDoubleArray2D((sim->levelOneChannels)[s], "Level 1 channel", qcode->nlogs, qcode->nlogs);
	}
	PrintDoubleArray2D((sim->levelOneCosets), "Level 1 coset probabilities", qcode->nstabs, qcode->nlogs);

	ConstructCumulative(sim->levelOneSynds, sim->levelOneCumul, qcode->nstabs);
	// Compute the importance distribution for level-1 if necessary.
	double *searchin = malloc(sizeof(double) * 2);
	if (sim->importance == 1) {
		searchin[0] = 0;
		searchin[1] = 1;
		double expo = PowerSearch(sim->syndprobs, qcode->nstabs, sim->outlierprobs, searchin);
		ConstructImportanceDistribution(sim->syndprobs, sim->levelOneImpDist, qcode->nstabs, expo);
		ConstructCumulative(sim->levelOneImpDist, sim->levelOneImpCumul, qcode->nstabs);
	}
	free(searchin);
	// Running averages have no meaning since we compute the exact average for
	// level 1. However, we'll set all the running averages to the exact one in
	// this case. int m; for (m = 0; m < sim->nmetrics; m ++){
	// (sim->runavg)[m][0] = (double) (sim->nbreaks); 	for (i = 0; i <
	// sim->nbreaks; i ++){ 		(sim->runstats)[i] = 1;
	// (sim->runavg)[m][i
	// + 1] = (sim->metricValues)[1][m];
	// 	}
	// }

	// copy lookup from sim to qcode before freeing
	if(copy_lookup==1)
	{
		for (s = 0; s < qcode->nstabs; s++)
			(qcode->dclookup)[s]= (sim->corrections)[s];
	}

	FreeSimParamsQECC(sim, qcode->N, qcode->K);
}

void ComputeLogicalChannels(struct simul_t **sims, struct qecc_t **qcode, struct constants_t *consts, double *****channels) {
	// Compute a logical channel for the required concatenation level.
	// The logical channel at a concatenated level l depends on N channels from
	// the previous concatenation level, and so on... until 7^l physical channels.
	// printf("Function: ComputeLogicalChannels\n");
	int *nphys = malloc(sizeof(int) * sims[0]->nlevels);
	int l;
	for (l = 0; l < sims[0]->nlevels; l++)
		nphys[l] = qcode[l]->N;
	// PrintIntArray1D(nphys, "nphys", sims[0]->nlevels);
	int *chans = malloc(sizeof(int) * (sims[0]->nlevels + 1));
	CountIndepLogicalChannels(chans, nphys, sims[0]->nlevels);
	// PrintIntArray1D(chans, "chans", sims[0]->nlevels + 1);

	// At every level, select a set of n channels, consider them as physical
	// channels and perform qcode to output a logical channel. Place this logical
	// channel in the channels array, at the succeeding level. To start with, we
	// will only initialize the last level with samples of the level-1 channels.
	int s, b, i, j, q, randsynd;
	double bias = 1, history = 1, expo;
	double *searchin = malloc(sizeof(double) * 2);
	double *impdist = malloc(sizeof(double) * qcode[0]->nstabs);
	double *impcumul = malloc(sizeof(double) * qcode[0]->nstabs);
	double ****inputchannels = malloc(sizeof(double ***));
	int *isPauli = malloc(sizeof(int) * (int)(1 + (int)(sims[0]->decoders[0] == 2)));

	// printf("Computing logical channels for %d levels.\n", sims[0]->nlevels);

	for (l = 1; l < sims[0]->nlevels; l++){
		// Allocate memory for the simulation parameters which depend on the error correcting code
		for (s = 0; s < 1 + (int)((sims[0]->decoders)[0] == 2); s++)
			AllocSimParamsQECC(sims[s], qcode[l]->N, qcode[l]->K);

		// Allocate memory for inputchannels
		inputchannels = (double ****)realloc(inputchannels, sizeof(double ***) * qcode[l]->N);
		MemManageInputChannels(inputchannels, qcode[l]->N, qcode[l]->nlogs, (sims[0]->decoders)[0], 0);

		// Perform coarsegraining of logical channels
		Coarsegrain(l - 1, sims, channels, chans[l], qcode[l]->nlogs);

		for (b = 0; b < chans[l + 1]; b++) {
		// printf("batch = %d of %d\n", b, chans[l + 1]);
		// Load the input channels on to the simulation structures and perform QECC.
		// Iterating in reverse order to accommodate partial ML decoder
		for (s = (int)((sims[0]->decoders)[0] == 2); s>=0; s--) {
			bias = 1;
			history = 1;
			isPauli[s] = 1;
			for (q = 0; q < qcode[l]->N; q++) {
				// printf("Simulation channel number %d, qubit %d:\n",s, q);
				for (i = 0; i < qcode[l]->nlogs; i++){
					for (j = 0; j < qcode[l]->nlogs; j++)
						(sims[s]->virtchan)[q][i][j] = channels[l - 1][qcode[l]->N * b + q][s][i][j];
					(sims[s]->pauli_probs)[q][i] = channels[l - 1][qcode[l]->N * b + q][s][qcode[l]->nlogs + 1][i];
				}

				if (isPauli[s] > 0)
					isPauli[s] = isPauli[s] * IsDiagonal((sims[s]->virtchan)[q], qcode[l]->nlogs);

				// PrintDoubleArray2D((sims[s]->virtchan)[q], "virtual channel", qcode[l]->nlogs, qcode[l]->nlogs);
				// PrintDoubleArray1D((sims[s]->pauli_probs)[q], "coset probabilities of chosen channel", qcode[l]->nlogs);
				bias *= channels[l - 1][qcode[l]->N * b + q][s][qcode[l]->nlogs][0];
				history *= channels[l - 1][qcode[l]->N * b + q][s][qcode[l]->nlogs][1];
				// printf("Syndrome probability: %g, bias = %g.\n", channels[l - 1][qcode[l]->N * b + q][s][qcode[l]->nlogs][1], channels[l - 1][qcode[l]->N * b + q][s][qcode[l]->nlogs][0]);
				// printf("======\n");
			}
			// printf("Going to perform SingleShotErrorCorrection on s = %d, isPauli = %d and frame = %d.\n", s, isPauli[s], (sims[s]->frames)[l]);
			// Pass minimum weight for main channel as decoding algo if partial ML decoding on
			if ((s == 0) && (sims[0]->decoders)[0] == 2)
				SingleShotErrorCorrection(isPauli[s], 0, 1, (sims[s]->frames)[l], qcode[l], sims[s], consts, 0);
			else if (s == 1)
				SingleShotErrorCorrection(isPauli[s], 0, 0, (sims[s]->frames)[l], qcode[l], sims[s], consts, 0);
			else
				SingleShotErrorCorrection(isPauli[s], 0, (sims[s]->decoders)[0], (sims[s]->frames)[l], qcode[l], sims[s], consts, 0);

			// If doing partial ML decoder for main channel,
			// copy the lookup table for the next simulation.
			if (s==1 && (sims[0]->decoders)[0] == 2){
				for (int synd = 0; synd < qcode[l]->nstabs; synd++){
					// printf("Copying lookup inside compute logical for synd %d, correction = %d \n",synd,(sims[1]->corrections)[synd]);
					(qcode[l]->dclookup)[synd]= (sims[1]->corrections)[synd];
				}
			}
			// Reference channel does ML while the original channel does lookup with the table generated by ML on the reference.
		}

		UpdateMetrics(l, bias, history, 0, qcode[l], sims[0], consts);

		if (l < (sims[0]->nlevels - 1)) {
			// Implementing partial ML decoder case
			// Main channel draws syndrome and passes to auxillary
			if ((sims[0]->decoders)[0] == 2){

				if (sims[0]->importance == 1){
					searchin[0] = 0;
					searchin[1] = 1;
					expo = PowerSearch(sims[0]->syndprobs, qcode[l]->nstabs, sims[0]->outlierprobs, searchin);
					ConstructImportanceDistribution(sims[0]->syndprobs, impdist, qcode[l]->nstabs, expo);
					ConstructCumulative(impdist, impcumul, qcode[l]->nstabs);
					randsynd = SampleCumulative(impcumul, qcode[l]->nstabs);
				}
				else
					randsynd = SampleCumulative(sims[0]->cumulative, qcode[l]->nstabs);
				// printf("Random syndrome = %d\n", randsynd);
				for (i = 0; i < qcode[l]->nlogs; i++)
					for (j = 0; j < qcode[l]->nlogs; j++)
					  for (s = 0; s < 2; s++)
							channels[l][b][s][i][j] = (sims[s]->effprocess)[randsynd][i][j];
				if(sims[0]->importance == 1){
					channels[l][b][0][qcode[l]->nlogs][0] = (sims[0]->syndprobs)[randsynd] / impdist[randsynd];
					channels[l][b][0][qcode[l]->nlogs][1] = history *(sims[0]->syndprobs)[randsynd];
					channels[l][b][0][qcode[l]->nlogs][2] = (sims[0]->syndprobs)[randsynd];

					channels[l][b][1][qcode[l]->nlogs][0] = (sims[1]->syndprobs)[randsynd] / impdist[randsynd];
					channels[l][b][1][qcode[l]->nlogs][1] = history *(sims[1]->syndprobs)[randsynd];
					channels[l][b][1][qcode[l]->nlogs][2] = (sims[1]->syndprobs)[randsynd];
				}
				else{
				channels[l][b][0][qcode[l]->nlogs][0] = 1;
				channels[l][b][0][qcode[l]->nlogs][1] = sims[0]->syndprobs[randsynd];
				channels[l][b][0][qcode[l]->nlogs][2] = sims[0]->syndprobs[randsynd];
				// 2. Drawing syndromes for the auxillary channel according to the
				// noisy channel syndrome distribution.
				channels[l][b][1][qcode[l]->nlogs][0] = 1;
				channels[l][b][1][qcode[l]->nlogs][1] = sims[1]->syndprobs[randsynd];
				channels[l][b][1][qcode[l]->nlogs][2] = sims[1]->syndprobs[randsynd];
				}
			}
			else if (sims[0]->importance == 0) {
				randsynd = SampleCumulative(sims[0]->cumulative, qcode[l]->nstabs);
				// printf("Random syndrome = %d\n", randsynd);
				for (i = 0; i < qcode[l]->nlogs; i++)
					for (j = 0; j < qcode[l]->nlogs; j++)
						channels[l][b][0][i][j] = (sims[0]->effprocess)[randsynd][i][j];
				channels[l][b][0][qcode[l]->nlogs][0] = 1;
				channels[l][b][0][qcode[l]->nlogs][1] = history * sims[0]->syndprobs[randsynd];
				channels[l][b][0][qcode[l]->nlogs][2] = sims[0]->syndprobs[randsynd];
			}
			else if (sims[0]->importance == 1) {
				// Compute a probability distribution where the probability of every
				// syndrome is given by a power of the original syndrome distribution.
				// The new distribution Q(s) is given by Eq. 6 of the article.
				// Sample a syndrome according to Q(s) and add a bias P(s)/Q(s).
				searchin[0] = 0;
				searchin[1] = 1;
				expo = PowerSearch(sims[0]->syndprobs, qcode[l]->nstabs, sims[0]->outlierprobs, searchin);
				ConstructImportanceDistribution(sims[0]->syndprobs, impdist, qcode[l]->nstabs, expo);
				ConstructCumulative(impdist, impcumul, qcode[l]->nstabs);
				randsynd = SampleCumulative(impcumul, qcode[l]->nstabs);
				// printf("randsynd = %d\n", randsynd);
				for (i = 0; i < qcode[l]->nlogs; i++)
					for (j = 0; j < qcode[l]->nlogs; j++)
						channels[l][b][0][i][j] = (sims[0]->effprocess)[randsynd][i][j];
				// printf("Populated channels.\n");
				channels[l][b][0][qcode[l]->nlogs][0] = (sims[0]->syndprobs)[randsynd] / impdist[randsynd];
				channels[l][b][0][qcode[l]->nlogs][1] = history * (sims[0]->syndprobs)[randsynd];
				channels[l][b][0][qcode[l]->nlogs][2] = (sims[0]->syndprobs)[randsynd];
			}
			else
				continue;
			}
		}
		// Free memory for inputchannels
		MemManageInputChannels(inputchannels, qcode[l]->N, qcode[l]->nlogs, (sims[0]->decoders)[0], 1);
		// Free simulation parameters that depend on the qcode
		for (s = 0; s < 1 + (int)(sims[0]->decoders[l] == 2); s++)
			FreeSimParamsQECC(sims[s], qcode[l]->N, qcode[l]->K);
	}
	// printf("Freeing all memory.\n");
	// Free memory
	free(searchin);
	free(isPauli);
	free(impdist);
	free(impcumul);
	free(chans);
	free(nphys);
	// printf("Freeing inputchannels.\n");
	// for (q = 0; q < qcode[sims[0]->nlevels - 1]->N; q ++)
	// 	free(inputchannels[q]);
	free(inputchannels);
}

void Performance(struct qecc_t **qcode, struct simul_t **sims, struct constants_t *consts) {
	// Compute logical error rates for a concatenation level.
	init_genrand(time(NULL)); // See the random number generator
	/*
	channels = list of l arrays.
	channels[i] = list of list of 7^(l-i-1) arrays.
	channels[i][b] = list of 3 arrays -- each
	corresponding to one type of sampling method. If importance sampling is
	turned off, there is only one array in the list. channels[i][b][s] = 4x4
	matrix denoting a level-i channel.
	*/
	int *nphys = malloc(sims[0]->nlevels * sizeof(int));
	int *nencs = malloc(sims[0]->nlevels * sizeof(int));
	int l;
	for (l = 0; l < sims[0]->nlevels; l++) {
		nphys[l] = qcode[l]->N;
		nencs[l] = qcode[l]->K;
	}
	// printf("Allocate memory to channels.\n");
	double *****channels = malloc(sizeof(double ****) * (sims[0]->nlevels));
	int nchans = MemManageChannels(channels, nphys, nencs, sims[0]->nlevels, (sims[0]->decoders)[0], 0);
	// Compute level-0 metrics.
	ComputeLevelZeroMetrics(sims[0], qcode[0]->N, qcode[0]->nlogs, consts);
	// Compute level-1 effective channels and syndromes.
	// PrintIntArray1D((sims[0]->decoders), "Performance: Decoders", sims[0]->nlevels);
	int s, decoder_to_use, copy_lookup;
	// Iterating reverse order to accommodate partial ML decoder
	for ( s = (int)((sims[0]->decoders)[0] == 2);s >= 0; s--)
	{
		copy_lookup=0; // whether to copy the lookup table at the end of Computing level one channels
		if(s==0 && (sims[0]->decoders)[0] == 2)
			decoder_to_use = 1;
		else if(s==1)
		{
			copy_lookup = 1;
		  decoder_to_use = 0;
		}
		else
			decoder_to_use = (sims[s]->decoders)[0];
		// printf("Computing level one channels in performance for channel %d\n", s);
		ComputeLevelOneChannels(sims[s], qcode[0], consts, decoder_to_use, copy_lookup);
	}
	// printf("Finished level-1 computations with %d channels.\n", (nchans));
	// PrintDoubleArray2D((sims[0]->logical)[1], "Logical channel", qcode[0]->nlogs, qcode[0]->nlogs);

	int c, i, j, m, randsynd;
	long t;
	if (sims[0]->nlevels > 1) {
		for (t = 0; t < sims[0]->nstats; t++) {
			// Fill the lowest level of the channels array with "nchans" samples of level-1 channels.
			// printf("Stat %ld, nchans = %d.\n", t, nchans);
			for (c = 0; c < nchans; c++) {
				if ((sims[0]->decoders)[0] == 2){
					if(sims[0]->importance == 1)
						randsynd = SampleCumulative(sims[0]->levelOneImpCumul, qcode[0]->nstabs);
					else
						randsynd = SampleCumulative(sims[0]->levelOneCumul, qcode[0]->nstabs);
					for (s = 0; s < 2; s++)
						for (i = 0; i < qcode[0]->nlogs; i++){
							for (j = 0; j < qcode[0]->nlogs; j++)
								channels[0][c][s][i][j] = (sims[s]->levelOneChannels)[randsynd][i][j];
							channels[0][c][s][1 + qcode[0]->nlogs][i] = (sims[s]->levelOneCosets)[randsynd][i];
						}
					if(sims[0]->importance == 1){
						channels[0][c][0][qcode[0]->nlogs][0] = (sims[0]->levelOneSynds)[randsynd] / (sims[0]->levelOneImpDist)[randsynd];
						channels[0][c][1][qcode[0]->nlogs][0] = (sims[1]->levelOneSynds)[randsynd] / (sims[1]->levelOneImpDist)[randsynd];
					}
					else{
						channels[0][c][0][qcode[0]->nlogs][0] = 1.0;
						channels[0][c][1][qcode[0]->nlogs][0] = 1.0;
					}

					channels[0][c][0][qcode[0]->nlogs][1] = (sims[0]->levelOneSynds)[randsynd];
					channels[0][c][0][qcode[0]->nlogs][2] = (sims[0]->levelOneSynds)[randsynd];

					channels[0][c][1][qcode[0]->nlogs][1] = (sims[1]->levelOneSynds)[randsynd];
					channels[0][c][1][qcode[0]->nlogs][2] = (sims[1]->levelOneSynds)[randsynd];
				}
				else if (sims[0]->importance == 0) {
					// Direct sampling
					randsynd = SampleCumulative(sims[0]->levelOneCumul, qcode[0]->nstabs);
					// printf("Random syndrome = %d, with probability: %g.\n", randsynd,
					// (sims[0]->levelOneSynds)[randsynd]);
					// printf("randsynd = %d\n", randsynd);
					for (i = 0; i < qcode[0]->nlogs; i++){
						for (j = 0; j < qcode[0]->nlogs; j++)
							channels[0][c][0][i][j] = (sims[0]->levelOneChannels)[randsynd][i][j];
						channels[0][c][0][1 + qcode[0]->nlogs][i] = (sims[0]->levelOneCosets)[randsynd][i];
					}
					channels[0][c][0][qcode[0]->nlogs][0] = 1.0;
					channels[0][c][0][qcode[0]->nlogs][1] = (sims[0]->levelOneSynds)[randsynd];
					channels[0][c][0][qcode[0]->nlogs][2] = (sims[0]->levelOneSynds)[randsynd];
				}
				else if (sims[0]->importance == 1) {
					// Draw a syndrome from the importance distribution specified by the
					// power-law scaling.
					randsynd = SampleCumulative(sims[0]->levelOneImpCumul, qcode[0]->nstabs);
					// printf("randsynd = %d\n", randsynd);
					for (i = 0; i < qcode[0]->nlogs; i++){
						for (j = 0; j < qcode[0]->nlogs; j++)
							channels[0][c][0][i][j] = (sims[0]->levelOneChannels)[randsynd][i][j];
						channels[0][c][0][1 + qcode[0]->nlogs][i] = (sims[0]->levelOneCosets)[randsynd][i];
					}
					channels[0][c][0][qcode[0]->nlogs][0] = (sims[0]->levelOneSynds)[randsynd] / (sims[0]->levelOneImpDist)[randsynd];
					channels[0][c][0][qcode[0]->nlogs][1] = (sims[0]->levelOneSynds)[randsynd];
					channels[0][c][0][qcode[0]->nlogs][2] = (sims[0]->levelOneSynds)[randsynd];
				}
				else
					continue;
			}

			// Compute average logical channels and average logical error rates.
			ComputeLogicalChannels(sims, qcode, consts, channels);
			for (s = 0; s < 1 + (int)((sims[0]->decoders)[0] == 2); s++){
				if (IsElement(sims[s]->runstats, sims[s]->nbreaks, t + 1) == 1) {
					for (m = 0; m < sims[0]->nmetrics; m++) {
						(sims[s]->runavg)[m][0] += 1;
						(sims[s]->runavg)[m][(int)(sims[s]->runavg)[m][0]] = (sims[s]->metricValues)[sims[s]->nlevels][m] / ((double)(t + 1));
					}
				}
			}
		}
	}

	// Normalize the average metrics.
	for (l = 1; l < sims[0]->nlevels; l++)
		UpdateMetrics(l, 1.0, 1.0, 1, qcode[l], sims[0], consts);

	// printf("Updated metrics\n");
	// PrintDoubleArray2D((sims[0]->logical)[1], "logical level-1 channel",
	// qcode[0]->nlogs, qcode[0]->nlogs);
	// PrintDoubleArray1D((sims[0]->metricValues)[1], "logical level-1 metrics",
	// sims[0]->nmetrics); PrintDoubleArray2D((sims[0]->logical)[1], "Logical
	// channel", qcode[0]->nlogs, qcode[0]->nlogs);

	// Free memory for channels.
	MemManageChannels(channels, nphys, nencs, sims[0]->nlevels, (sims[0]->decoders)[0], 1);
	free(channels);
	free(nphys);
	free(nencs);

	// printf("Freed channels\n");
	// printf("Done Performance.\n");
}
