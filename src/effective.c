#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "../mt19937/mt19937ar.h"
#include "constants.h"
#include "printfuns.h"
#include "utils.h"
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

void UpdateMetrics(int level, long double bias, long double history, int isfinal, struct qecc_t *qcode, struct simul_t *sim, struct constants_t *consts) {
	// Compute metrics for all the effective channels and update the average value
	// of the metrics. Metric values.
	// printf("Updating metrics.\n");
	double normalization = 0;
	double *metvals = malloc(sim->nmetrics * sizeof(double));
	long double *avg = malloc((sim->nmetrics + qcode->nlogs * qcode->nlogs) * sizeof(long double));
	int m;
	for (m = 0; m < sim->nmetrics + qcode->nlogs * qcode->nlogs; m++)
		avg[m] = 0;

	int s, i, j;
	if (isfinal == 0) {
		for (s = 0; s < qcode->nstabs; s++) {
			if ((sim->syndprobs)[s] > consts->atol) {
				// Compute metrics.
				ComputeMetrics(metvals, sim->nmetrics, sim->metricsToCompute, (sim->effprocess)[s], sim->chname, consts);
				for (m = 0; m < sim->nmetrics; m++)
					avg[m] += (long double) metvals[m] * (sim->syndprobs)[s];
				// Compute average channel.
				for (i = 0; i < qcode->nlogs; i++)
					for (j = 0; j < qcode->nlogs; j++)
						avg[sim->nmetrics + i * qcode->nlogs + j] += (sim->effprocess)[s][i][j] * (sim->syndprobs)[s];
			}
		}
		/*
		if (level == 2){
			printf("Level = %d, Bias = %.3Le\n", level, bias);
			PrintLongDoubleArray1D(avg, "Avg Metric values and channels", sim->nmetrics + qcode->nlogs * qcode->nlogs);
		}
		*/
		
		// printf("Level = %d, Bias = %.3Le\n", level, bias);
		// PrintLongDoubleArray1D(avg, "Avg Metric values and channels", sim->nmetrics + qcode->nlogs * qcode->nlogs);

		// Average of metrics.
		for (m = 0; m < sim->nmetrics; m++) {
			(sim->metricValues)[level + 1][m] += bias * avg[m];
			(sim->sumsq)[level + 1][m] += (double) powl(bias * avg[m], 2);
		}
		// Populating average logical channel.
		for (i = 0; i < qcode->nlogs; i++) {
			for (j = 0; j < qcode->nlogs; j++) {
				(sim->logical)[level + 1][i][j] += bias * (double) avg[sim->nmetrics + i * qcode->nlogs + j];
				(sim->sumsq)[level + 1][sim->nmetrics + i * qcode->nlogs + j] += (double) powl(bias * avg[sim->nmetrics + i * qcode->nlogs + j], 2);
			}
		}
		
		// Syndrome-metric binning.
		int mbin, sbin = GetBinPosition(fabsl(history), sim->nbins, sim->maxbin);
		for (m = 0; m < sim->nmetrics; m++) {
			mbin = GetBinPosition(fabsl(avg[m]), sim->nbins, sim->maxbin);
			(sim->bins)[level + 1][m][sbin][mbin] ++;
		}
		
		// Update the number of statistics done for the level.
		// printf("totalbias[%d] = %Lf\n", level + 1, (sim->totalbias)[level + 1]);
		(sim->statsperlevel)[level + 1] ++;
		// Force normalization of level - 3 simulations
		// (sim->totalbias)[level + 1] += bias;
	}
	else{
		// After all the simulations are done, the average metrics are to be computed by diving the metricValues by the total number of statistics done for that level.
		// Note that sampling the level - 2 syndromes (and also for higher levels) does not have access to the full distribution (or, even approximately the full distribution, for low sample sizes).
		// So we need to apply a self-normalized version of importance sampling wherein we normalize the estimator by the sum of all bias values.
		/*
		if (level < 2)
			normalization = (double) ((sim->statsperlevel)[level + 1]);
		else
			normalization = (double) (sim->totalbias)[level + 1];
		*/
		normalization = (double) ((sim->statsperlevel)[level + 1]);

		// printf("Level %d\nTotal infidelity = %.8e and normalization = %.4e.\n", level, normalization, (sim->metricValues)[level + 1][0]);

		// Normalize metrics
		for (m = 0; m < sim->nmetrics; m++) {
			(sim->metricValues)[level + 1][m] /= normalization;
			(sim->variance)[level + 1][m] = 1 / (normalization * (normalization - 1)) * (sim->sumsq)[level + 1][m] - pow((sim->metricValues)[level + 1][m], 2);
		}
		// Normalize logical channel
		for (i = 0; i < qcode->nlogs; i++) {
			for (j = 0; j < qcode->nlogs; j++) {
				(sim->logical)[level + 1][i][j] /= normalization;
				(sim->variance)[level + 1][sim->nmetrics + i * qcode->nlogs + j] = 1 / (normalization * (normalization - 1)) * sim->sumsq[level + 1][sim->nmetrics + i * qcode->nlogs + j] - pow((sim->logical)[level + 1][i][j], 2);
			}
		}
		if (_IsChannel((sim->logical)[level + 1], consts, consts->atol, 1, 0) == 0){
			printf("!!!!!!!!!!\n");
			printf("Invalid average logical channel at level %d.\n", level + 1);
			printf("!!!!!!!!!!\n");
		}
	}
	// Free memory.
	free(metvals);
	free(avg);
}


void ComputeLevelOneChannels(struct simul_t *sim, struct qecc_t *qcode, struct constants_t *consts, int decoder, int copy_lookup) {
	/*
		Compute the effective channels and syndrome probabilities for all level-1 syndromes.
		The pre-computation is to avoid re-computing level-1 syndromes for every new top-level syndrome.
		Load the physical channels on to the simulation structure and perform single shot quantum error correction.
	*/
	// printf("Function: ComputeLevelOneChannels with decoder %d\n", decoder);
	int q, i, j, isPauli = 1, synd_threshold;
	synd_threshold = Min(OrderOfMagnitude(consts->min_syndprob, 10), OrderOfMagnitude(1/((long double) sim->nstats), 10));
	AllocSimParamsQECC(sim, qcode->N, qcode->K);

	// For the topmost level, we don't need to compute conditional channels accurately.
	// Here, we will compute E_s * P_s as a single entity, i.e., we will not normalize E_s * P_s by P_s like we do usually for the lower levels.
	if (sim->nlevels == 1)
		sim->skipsyndromes = 1;

	if ((sim->decoders)[0] == 3){
		ComputeCosetProbsLevelOne(sim->mpinfo, qcode->nlogs, qcode->nstabs, sim->cosetprobs);
		// PrintLongDoubleArray2D(sim->cosetprobs, "Coset probabilities", 64, 4);
	}

	if ((sim->iscorr == 0) || (sim->iscorr == 2)) {
		for (q = 0; q < qcode->N; q++) {
			for (i = 0; i < qcode->nlogs; i++){
				for (j = 0; j < qcode->nlogs; j++){
					if (sim->iscorr == 0)
						(sim->virtchan)[q][i][j] = (long double) ((sim->physical)[i * qcode->nlogs + j]);
					else
						(sim->virtchan)[q][i][j] = (long double) ((sim->physical)[q * qcode->nlogs * qcode->nlogs + i * qcode->nlogs + j]);
				}
			}
			if (isPauli > 0)
				isPauli = isPauli * IsDiagonal((sim->virtchan)[q], qcode->nlogs);
			// printf("Qubit %d.\n", q);
			// PrintDoubleArrayDiag((sim->virtchan)[q], "Level 0 channel", qcode->nlogs);
		}

		// printf("Simulating an uncorrelated channel: iscorr = %d, skipsyndromes = %d, synd_threshold = %d and decoder %d.\n", sim->iscorr, sim->skipsyndromes, synd_threshold, decoder);
		SingleShotErrorCorrection(isPauli, sim->iscorr, decoder, (sim->frames)[0], qcode, sim, consts, 1, synd_threshold);
	}
	else{
		if (sim->iscorr == 1)
			isPauli = 1;
		else
			isPauli = 0;
		// This runs for iscorr = 1 as well as iscorr = 3.
		// Correlated channel -- the physical channel contains the full process matrix
		// printf("Simulating a correlated channel: level: 1, isPauli = %d, iscorr = %d, skipsyndromes = %d, synd_threshold = %d and decoder %d.\n", isPauli, sim->iscorr, sim->skipsyndromes, synd_threshold, decoder);
		SetFullProcessMatrix(qcode, sim, sim->physical, isPauli);
		// printf("Running SingleShotErrorCorrection.\n");
		SingleShotErrorCorrection(isPauli, sim->iscorr, decoder, (sim->frames)[0], qcode, sim, consts, 1, synd_threshold);
	}
	// printf("Completed SingleShotErrorCorrection for level 1.\n");

	UpdateMetrics(0, 1, 1, 0, qcode, sim, consts);
	int s;

	for (s = 0; s < qcode->nstabs; s++) {
		(sim->levelOneSynds)[s] = (sim->syndprobs)[s];
		for (i = 0; i < qcode->nlogs; i++){
			for (j = 0; j < qcode->nlogs; j++)
				(sim->levelOneChannels)[s][i][j] = (sim->effprocess)[s][i][j];
			(sim->levelOneCosets)[s][i] = (sim->cosetprobs)[s][i];
		}
		// printf("s = %d: P(s) = %.15Lf,\n", s, (sim->levelOneSynds)[s]);
		// PrintLongDoubleArray1D((sim->levelOneCosets)[s], "Level 1 coset probability", qcode->nlogs);
		// PrintLongDoubleArray2D((sim->levelOneChannels)[s], "Level 1 channel", qcode->nlogs, qcode->nlogs);
	}
	// PrintDoubleArray2D((sim->levelOneCosets), "Level 1 coset probabilities", qcode->nstabs, qcode->nlogs);
	ConstructCumulative(sim->levelOneSynds, sim->levelOneCumul, qcode->nstabs);
	// Compute the importance distribution for level-1 if necessary.
	double *searchin = malloc(sizeof(double) * 2);
	if (sim->importance == 1) {
		searchin[0] = 0;
		searchin[1] = 1;
		// Method 1:
		SetOutlierProbs(sim->infidelity, qcode->D, 1, sim->outlierprobs);
		double expo = PowerSearch(sim->syndprobs, qcode->nstabs, sim->outlierprobs, searchin);
		// Method 2:
		// double expo = SetExponent(sim->infidelity, qcode->D, qcode->N, 2, sim->levelOneSynds);
		// printf("level-1 exponent = %g.\n", expo);
		ConstructImportanceDistribution(sim->syndprobs, sim->levelOneImpDist, qcode->nstabs, expo);
		ConstructCumulative(sim->levelOneImpDist, sim->levelOneImpCumul, qcode->nstabs);
	}
	free(searchin);
	// Running averages have no meaning since we compute the exact average for level 1.

	// copy lookup from sim to qcode before freeing
	if(copy_lookup==1)
	{
		for (s = 0; s < qcode->nstabs; s++)
			(qcode->dclookup)[s]= (sim->corrections)[s];
	}

	FreeSimParamsQECC(sim, qcode->N, qcode->K);
}


void ComputeLogicalChannels(struct simul_t **sims, struct qecc_t **qcode, struct constants_t *consts, long double *****channels) {
	/*
		Compute a logical channel for the required concatenation level.
		The logical channel at a concatenated level l depends on N channels from
		the previous concatenation level, and so on... until N^l physical channels.
	*/
	int *nphys = malloc(sizeof(int) * sims[0]->nlevels);
	int l;
	for (l = 0; l < sims[0]->nlevels; l++)
		nphys[l] = qcode[l]->N;
	
	int *chans = malloc(sizeof(int) * (sims[0]->nlevels + 1));
	CountIndepLogicalChannels(chans, nphys, sims[0]->nlevels);
	
	/*
		At every level, select a set of n channels, consider them as physical channels and perform qcode to output a logical channel.
		Place this logical channel in the channels array, at the succeeding level.
		To start with, we will only initialize the last level with samples of the level-1 channels.
	*/
	int s, b, i, j, q, randsynd = 0, history_order, synd_threshold, threshold;
	long double bias = 1, history = 1;
	double expo;
	double *searchin = malloc(sizeof(double) * 2);
	long double *impdist = malloc(sizeof(long double) * qcode[0]->nstabs);
	long double *impcumul = malloc(sizeof(long double) * qcode[0]->nstabs);
	long double ****inputchannels = malloc(sizeof(long double ***));
	int *isPauli = malloc(sizeof(int) * (int)(1 + (int)(sims[0]->decoders[0] == 2)));

	for (l = 1; l < sims[0]->nlevels; l++){
		// Allocate memory for the simulation parameters which depend on the error correcting code
		for (s = 0; s < 1 + (int)((sims[0]->decoders)[0] == 2); s++)
			AllocSimParamsQECC(sims[s], qcode[l]->N, qcode[l]->K);

		// For the topmost level, we don't need to compute conditional channels accurately.
		// Here, we will compute E_s * P_s as a single entity, i.e., we will not normalize E_s * P_s by P_s like we do usually for the lower levels.
		if (l == sims[0]->nlevels - 1)
			sims[0]->skipsyndromes = 1;

		// Allocate memory for inputchannels
		inputchannels = (long double ****)realloc(inputchannels, sizeof(long double ***) * qcode[l]->N);
		MemManageInputChannels(inputchannels, qcode[l]->N, qcode[l]->nlogs, (sims[0]->decoders)[0], 0);

		// Perform coarsegraining of logical channels
		// Coarsegrain(l - 1, sims, channels, chans[l], qcode[l]->nlogs);

		for (b = 0; b < chans[l]; b++) {
			// Load the input channels on to the simulation structures and perform QECC.
			// Iterating in reverse order to accommodate partial ML decoder
			for (s = (int)((sims[0]->decoders)[0] == 2); s>=0; s--) {
				bias = 1;
				history = 1;
				isPauli[s] = 1;
				// printf("Simulating QECC for a code with N = %d, nstabs = %d\n", qcode[l]->N, qcode[l]->nstabs);
				for (q = 0; q < qcode[l]->N; q++) {
					for (i = 0; i < qcode[l]->nlogs; i++){
						for (j = 0; j < qcode[l]->nlogs; j++)
							(sims[s]->virtchan)[q][i][j] = channels[l - 1][qcode[l]->N * b + q][s][i][j];
						(sims[s]->pauli_probs)[q][i] = channels[l - 1][qcode[l]->N * b + q][s][qcode[l]->nlogs + 1][i];
					}
					// Check for Pauli channels in the simulation.
					if (isPauli[s] > 0)
						isPauli[s] = isPauli[s] * IsDiagonal((sims[s]->virtchan)[q], qcode[l]->nlogs);

					bias *= channels[l - 1][qcode[l]->N * b + q][s][qcode[l]->nlogs][0];
					history *= channels[l - 1][qcode[l]->N * b + q][s][qcode[l]->nlogs][1];
					
					// printf("Channel on Qubit %d\n", q);
					// PrintLongDoubleArray2D((sims[s]->virtchan)[q], "Channel", qcode[l]->nlogs, qcode[l]->nlogs);
				}
				// Compute the minimum syndrome probability that should be sampled.
				history_order = OrderOfMagnitude(history, 10);
				synd_threshold = OrderOfMagnitude(consts->min_syndprob, 10) - history_order + OrderOfMagnitude(1/((long double) sims[0]->nstats), 10);
				threshold = Min(Max((double) OrderOfMagnitude(consts->min_syndprob, 10), (double) synd_threshold), -5);
				
				// printf("Going to perform SingleShotErrorCorrection with decoder = %d, isPauli = %d, bias = %.5Le, skipsyndromes = %d, history ~ 10^%d, and synd_threshold = %d.\n", (sims[0]->decoders)[0], isPauli[s], bias, sims[s]->skipsyndromes, history_order, threshold);
				
				// Pass minimum weight for main channel as decoding algo if partial ML decoding on
				if (s == 0){
					if ((sims[0]->decoders)[0] == 2)
						SingleShotErrorCorrection(isPauli[s], 0, 1, (sims[s]->frames)[l], qcode[l], sims[s], consts, 0, threshold);
					else if (((sims[0]->decoders)[0] == 3) && (l < (sims[0]->nlevels - 1))) // Soft decoding -- should not apply logical operations at intermediate levels.
						SingleShotErrorCorrection(isPauli[s], 0, 4, (sims[s]->frames)[l], qcode[l], sims[s], consts, 0, threshold);
					else
						SingleShotErrorCorrection(isPauli[s], 0, (sims[s]->decoders)[0], (sims[s]->frames)[l], qcode[l], sims[s], consts, 0, threshold);
				}
				else if (s == 1)
					SingleShotErrorCorrection(isPauli[s], 0, 0, (sims[s]->frames)[l], qcode[l], sims[s], consts, 0, threshold);
				else{};

				// printf("Done SingleShotErrorCorrection with isPauli = %d, bias = %.5Le, skipsyndromes = %d, history ~ 10^%d, and synd_threshold = %d.\n", isPauli[s], bias, sims[s]->skipsyndromes, history_order, threshold);
				
				// If doing partial ML decoder for main channel, copy the lookup table for the next simulation.
				if ((s == 1) && ((sims[0]->decoders)[0] == 2))
					for (int synd = 0; synd < qcode[l]->nstabs; synd++)
						(qcode[l]->dclookup)[synd]= (sims[1]->corrections)[synd];
				// Reference channel does ML while the original channel does lookup with the table generated by ML on the reference.
			}

			// PrintLongDoubleArray1D(sims[s]->syndprobs, "Syndrome probabilities", qcode[l]->nstabs);

			UpdateMetrics(l, bias, history, 0, qcode[l], sims[0], consts);

			if (l < (sims[0]->nlevels - 1)) {
				if (sims[0]->importance == 0){
					randsynd = SampleCumulative(sims[0]->cumulative, qcode[l]->nstabs);
					bias *= 1;
					// impdist[randsynd] = (sims[0]->syndprobs)[randsynd];
					// printf("Done Sampling cumulative with randsynd = %d\n", randsynd);
					history *= (sims[0]->syndprobs)[randsynd];
					// printf("Done sampling\n");
				}
				else if (sims[0]->importance == 1){
					searchin[0] = 0;
					searchin[1] = 1;
					// Method 1:
					SetOutlierProbs(sims[0]->infidelity, qcode[0]->D, l+1, sims[0]->outlierprobs);
					expo = PowerSearch(sims[0]->syndprobs, qcode[l]->nstabs, sims[0]->outlierprobs, searchin);
					// Method 2:
					// expo = SetExponent(sims[0]->infidelity, qcode[0]->D, qcode[0]->N, l + 2, sims[0]->syndprobs);
					// printf("level %d exponent = %g.\n", l, expo);
					ConstructImportanceDistribution(sims[0]->syndprobs, impdist, qcode[l]->nstabs, expo);
					ConstructCumulative(impdist, impcumul, qcode[l]->nstabs);
					randsynd = SampleCumulative(impcumul, qcode[l]->nstabs);
					bias *= (sims[0]->syndprobs)[randsynd] / impdist[randsynd];
					history *= impdist[randsynd];
				}
				else{};
				// printf("Random syndrome: %d, from exponent: %g and bias = %Lf. (sims[0]->syndprobs)[randsynd] = %.15Lf, impdist[randsynd] = %.15Lf.\n", randsynd, expo, bias, (sims[0]->syndprobs)[randsynd], impdist[randsynd]);
				for (i = 0; i < qcode[l]->nlogs; i++){
					for (j = 0; j < qcode[l]->nlogs; j++)
						channels[l][b][0][i][j] = (sims[0]->effprocess)[randsynd][i][j];
					channels[l][b][0][1 + qcode[0]->nlogs][i] = (sims[0]->cosetprobs)[randsynd][i];
				}
				channels[l][b][0][qcode[l]->nlogs][0] = bias;
				channels[l][b][0][qcode[l]->nlogs][1] = history;
				channels[l][b][0][qcode[l]->nlogs][2] = (sims[0]->syndprobs)[randsynd];
			}
			// printf("Loaded channels for next run\n");
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
	free(inputchannels);
}

void Performance(struct qecc_t **qcode, struct simul_t **sims, struct constants_t *consts) {
	// Compute logical error rates for a concatenation level.
	init_genrand(time(NULL)); // See the random number generator
	int *nphys = malloc(sims[0]->nlevels * sizeof(int));
	int *nencs = malloc(sims[0]->nlevels * sizeof(int));
	int l;
	for (l = 0; l < sims[0]->nlevels; l++) {
		nphys[l] = qcode[l]->N;
		nencs[l] = qcode[l]->K;
	}

	/*
		channels = list of l arrays.
		channels[i] = list of list of 7^(l-i-1) arrays.
		channels[i][b] = list of 3 arrays -- each
		corresponding to one type of sampling method. If importance sampling is
		turned off, there is only one array in the list. channels[i][b][s] = 4x4
		matrix denoting a level-i channel.
	*/
	long double *****channels = malloc(sizeof(long double ****) * (sims[0]->nlevels));
	int nchans = MemManageChannels(channels, nphys, nencs, sims[0]->nlevels, (sims[0]->decoders)[0], 0);
	
	int s, decoder_to_use;
	int copy_lookup = 0; // whether to copy the lookup table at the end of Computing level one channels
	// Iterating reverse order to accommodate partial ML decoder
	for (s = (int)((sims[0]->decoders)[0] == 2); s >= 0; s--)
	{
		if (s == 0){
			if ((sims[0]->decoders)[0] == 2)
				decoder_to_use = 1;
			else if (((sims[0]->decoders)[0] == 3) && (sims[0]->nlevels > 1)) // Soft decoding -- should not apply logical operations at intermediate levels.
				decoder_to_use = 4;
			else
				decoder_to_use = (sims[s]->decoders)[0];
		}
		else if (s == 1){
			copy_lookup = 1;
			decoder_to_use = 0;
		}
		else{};
		ComputeLevelOneChannels(sims[s], qcode[0], consts, decoder_to_use, copy_lookup);
	}
	
	int c, i, j, m, randsynd;
	long double bias;
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
						
						channels[0][c][0][qcode[0]->nlogs][1] = (sims[0]->levelOneImpDist)[randsynd];
						channels[0][c][1][qcode[0]->nlogs][1] = (sims[1]->levelOneImpDist)[randsynd];
					}
					else{
						channels[0][c][0][qcode[0]->nlogs][0] = 1.0;
						channels[0][c][1][qcode[0]->nlogs][0] = 1.0;

						channels[0][c][0][qcode[0]->nlogs][1] = (sims[0]->levelOneSynds)[randsynd];
						channels[0][c][1][qcode[0]->nlogs][1] = (sims[1]->levelOneSynds)[randsynd];
					}
					channels[0][c][0][qcode[0]->nlogs][2] = (sims[0]->levelOneSynds)[randsynd];
					channels[0][c][1][qcode[0]->nlogs][2] = (sims[1]->levelOneSynds)[randsynd];
				}
				else{
					if (sims[0]->importance == 0){
						randsynd = SampleCumulative(sims[0]->levelOneCumul, qcode[0]->nstabs);
						bias = 1;
					}
					else{
						randsynd = SampleCumulative(sims[0]->levelOneImpCumul, qcode[0]->nstabs);
						bias = (sims[0]->levelOneSynds)[randsynd] / (sims[0]->levelOneImpDist)[randsynd];
					}
					// printf("Random syndrome: %d\n", randsynd);
					for (i = 0; i < qcode[0]->nlogs; i++){
						for (j = 0; j < qcode[0]->nlogs; j++)
							channels[0][c][0][i][j] = (sims[0]->levelOneChannels)[randsynd][i][j];
						channels[0][c][0][1 + qcode[0]->nlogs][i] = (sims[0]->levelOneCosets)[randsynd][i];
					}
					channels[0][c][0][qcode[0]->nlogs][0] = bias;
					if (sims[0]->importance == 0)
						channels[0][c][0][qcode[0]->nlogs][1] = (sims[0]->levelOneSynds)[randsynd];
					else
						channels[0][c][0][qcode[0]->nlogs][1] = (sims[0]->levelOneImpDist)[randsynd];
					channels[0][c][0][qcode[0]->nlogs][2] = (sims[0]->levelOneSynds)[randsynd];
				}
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

	// Free memory for channels.
	MemManageChannels(channels, nphys, nencs, sims[0]->nlevels, (sims[0]->decoders)[0], 1);
	free(channels);
	free(nphys);
	free(nencs);
	
	// printf("Done Performance.\n");
}
