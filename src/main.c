#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "mt19937/mt19937ar.h" // Random number generator
#include "rand.h"
#include "printfuns.h"
#include "utils.h"
#include "linalg.h"
#include "constants.h"
#include "memory.h"
#include "qecc.h"
#include "sampling.h"
#include "checks.h"
#include "logmetrics.h"
#include "benchmark.h"

int main(int argc, char **argv)
{
	/*
	This function is simply to test all the C functions in the converted/ folder.
	*/
	clock_t begin = clock();
	printf("bmark: source last compiled on %s at %s.\n", __DATE__, __TIME__);
	if (argc < 2){
		printf("Usage: ./bmark <file name> <function name>\n");
		return 0;
	}
	int i, j, k;
	double complex **mat = malloc(sizeof(double complex) * 4);
	for (i = 0; i < 4; i ++)
		mat[i] = malloc(sizeof(double complex) * 4);

	long double **matLD = malloc(sizeof(long double) * 4);
	for (i = 0; i < 4; i ++)
		matLD[i] = malloc(sizeof(long double) * 4);
	/* Creating a random number generator.
		See https://stackoverflow.com/questions/822323/how-to-generate-a-random-int-in-c on why not to use the in-built rand() function.
		Instead, we use the Mersenne Twister random number generator explained in http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/emt19937ar.html .
		The source code can be found in mt19937/ folder.
	*/
	init_genrand(time(NULL));
	// The file and the function to be tested are obtained from shell.
	char *file = malloc(sizeof(char) * 100);
	sprintf(file, "%s", argv[1]);
	char *func = malloc(sizeof(char) * 100);
	sprintf(func, "%s", argv[2]);

	if (strncmp(file, "sampling", 8) == 0){
		// Testing the functions in samplig.c.
		printf("File: samplig.c.\n");
		if (strncmp(func, "ConstructCumulative", 19) == 0){
			// Testing ConstructCumulative.
			printf("Testing ConstructCumulative.\n");
			int size = 10;
			// Creating a test cumulative distribution.
			long double *dist = malloc(sizeof(long double) * size);
			int i;
			long double sum = 0;
			for (i = 0; i < size; i ++){
				dist[i] = (long double) genrand_real3();
				sum += dist[i];
			}
			for (i = 0; i < size; i ++)
				dist[i] = dist[i]/sum;
			PrintLongDoubleArray1D(dist, "True Distribution", size);
			long double *cumul = malloc(sizeof(long double) * size);
			ConstructCumulative(dist, cumul, size);
			PrintLongDoubleArray1D(cumul, "Cumulative distribution", size);
			// free memory
			free(dist);
			free(cumul);
		}
		if (strncmp(func, "PowerSearch", 11) == 0){
			// Testing PowerSearch.
			printf("Function: PowerSearch.\n");
			int size = 10;
			// Creating a test cumulative distribution.
			long double *dist = malloc(sizeof(long double) * size);
			int i;
			dist[0] = (long double) genrand_real3();
			long double sum = dist[0];
			for (i = 1; i < size; i ++){
				dist[i] = 0.01 * (long double) genrand_real3();
				sum += dist[i];
			}
			for (i = 0; i < size; i ++)
				dist[i] = dist[i]/sum;
			double window[2] = {0.45, 0.5};
			double searchin[2] = {0, 1};
			PrintLongDoubleArray1D(dist, "Probability distribution: P", size);
			printf("Searching for k in [%g, %g] such that (P^k)[0] > %g.\n", searchin[0], searchin[1], window[0]);
			double exponent = PowerSearch(dist, size, window, searchin);
			printf("k = %g.\n", exponent);
			long double *powerdist = malloc(sizeof(long double) * size);
			ConstructImportanceDistribution(dist, powerdist, size, exponent);
			PrintLongDoubleArray1D(powerdist, "New distribution", size);
			// free memory
			free(dist);
			free(powerdist);
		}
	}

	if (strncmp(file, "constants", 9) == 0){
		// Testing the functions in constants.c.
		if (strncmp(func, "InitConstants", 13) == 0){
			// Testing the function InitConstants. This doesn't output anything.
			struct constants_t *consts = malloc(sizeof(struct constants_t));
			InitConstants(consts);
			PrintComplexArray2D(consts->choitochi, "Choi to Chi transformation matrix", 16, 16);
			FreeConstants(consts);
			free(consts);
		}
	}

	if (strncmp(file, "checks", 5) == 0){
		// Testing functions in checks.c.
		for (i = 0; i < 4; i ++)
			for (j = 0; j < 4; j ++)
				mat[i][j] = genrand_real3() + genrand_real3() * I;
		if (strncmp(argv[1], "IsPositive", 10) == 0){
			printf("Function: IsPositive.\n");
			PrintComplexArray2D(mat, "Matrix", 4, 4);
			int ispos = IsPositive(mat, 1E-12);
			if (ispos == 1)
				printf("is positive.\n");
			else{
				printf("is not positive, but\n");
				double complex **matP = malloc(sizeof(double complex) * 4);
				for (i = 0; i < 4; i ++){
					matP[i] = malloc(sizeof(double complex) * 4);
					for (j = 0; j < 4; j ++){
						matP[i][j] = 0;
						for (k = 0; k < 4; k ++)
							matP[i][j] += mat[i][k] * conj(mat[j][k]);
					}
				}
				PrintComplexArray2D(matP, "M . M^\\dag", 4, 4);
				ispos = IsPositive(matP, 1E-12);
				if (ispos == 1)
					printf("is positive.\n");
				else
					printf("Error in function.\n");
				// Free memory
				for (i = 0; i < 4; i ++)
					free(matP[i]);
				free(matP);
			}
		}
		if (strncmp(argv[1], "IsHermitian", 11) == 0){
			printf("Function: IsHermitian.\n");
			PrintComplexArray2D(mat, "Matrix M", 4, 4);
			int isherm = IsHermitian(mat, 1E-12);
			if (isherm == 1)
				printf("is Hermitian.\n");
			else{
				printf("is not Hermitian, but\n");
				double complex **matH = malloc(sizeof(double complex) * 4);
				for (i = 0; i < 4; i ++){
					matH[i] = malloc(sizeof(double complex) * 4);
					for (j = 0; j < 4; j ++)
						matH[i][j] = mat[i][j] + conj(mat[j][i]);
				}
				PrintComplexArray2D(matH, "M + M^\\dag", 4, 4);
				isherm = IsHermitian(matH, 1E-12);
				if (isherm == 1)
					printf("is Hermitian.\n");
				else
					printf("Error in function.\n");
				// Free memory
				for (i = 0; i < 4; i ++)
					free(matH[i]);
				free(matH);
			}
		}
		if (strncmp(argv[1], "IsTraceOne", 10) == 0){
			printf("Function: IsTraceOne.\n");
			PrintComplexArray2D(mat, "Matrix", 4, 4);
			int istrone = IsTraceOne(mat, 1E-12);
			if (istrone == 1)
				printf("has unit trace.\n");
			else{
				printf("does not have unit trace but\n");
				double trace = 0;
				for (i = 0; i < 4; i ++)
					trace = trace + creal(mat[i][i]);
				double complex **matN = malloc(sizeof(double complex) * 4);
				for (i = 0; i < 4; i ++){
					matN[i] = malloc(sizeof(double complex) * 4);
					for (j = 0; j < 4; j ++)
						matN[i][j] = mat[i][j]/trace;
					matN[i][i] = creal(matN[i][i]);
				}
				PrintComplexArray2D(matN, "M/trace", 4, 4);
				istrone = IsTraceOne(matN, 1E-12);
				if (istrone == 1)
					printf("has unit trace.\n");
				else
					printf("Error in function.\n");
				// Free memory
				for (i = 0; i < 4; i ++)
					free(matN[i]);
				free(matN);
			}
		}
		if (strncmp(argv[1], "IsState", 7) == 0){
			printf("Function: IsState.\n");
			PrintComplexArray2D(mat, "Matrix", 4, 4);
			int isstate = IsState(mat, 1E-12, 1, 0);
			if (isstate == 0)
				printf("is not a state.\n");
			else
				printf("is a state.\n");
		}
		if (strncmp(argv[1], "IsPDF", 5) == 0){
			long double sum = 0;
			long double *dist = malloc(sizeof(long double) * 100);
			for (i = 0; i < 100; i ++){
				dist[i] = (long double) genrand_real3();
				sum += dist[i];
			}
			PrintLongDoubleArray1D(dist, "Un-normalized distribution", 100);
			int ispdf = IsPDF(dist, 100);
			printf("has ispdf = %d.\n", ispdf);
			for (i = 0; i < 100; i ++)
				dist[i] = dist[i]/sum;
			PrintLongDoubleArray1D(dist, "And after normalization", 100);
			ispdf = IsPDF(dist, 100);
			printf("it has ispdf = %d.\n", ispdf);
			// Free memory
			free(dist);
		}
	}

	if (strncmp(file, "logmetrics", 10) == 0){
		// Testing functions in logmetrics.c.
		// REQUIRES SERIOUS REVISION DUE TO CHOI BEING REPLACED BY PTM.
		if (strncmp(func, "ComputeMetrics", 14) == 0){
			printf("Function: ComputeMetrics.\n");
			// Name of the channel
			char *chname = malloc(sizeof(char) * 100);
			sprintf(chname, "Random Channel");

			// Assign names of metrics whose values must be computed.
			int nmetrics = 6;
			char **metrics = malloc(sizeof(char *) * nmetrics);
			for (i = 0; i < nmetrics; i ++)
				metrics[i] = malloc(sizeof(char) * 100);
			sprintf(metrics[0], "frb");
			sprintf(metrics[1], "infid");
			sprintf(metrics[2], "processfidelity");
			sprintf(metrics[3], "trn");
			sprintf(metrics[4], "entropy");
			sprintf(metrics[5], "np1");

			// Create an error channel -- any complex matrix which is positive definte matrix and has unit trace, can be the input channel's Choi matrix.
			// Create a random matrix M and declare (M + M^dag)/(2 * trace(M)) to be the input Choi matrix.
			for (i = 0; i < 4; i ++)
				for (j = 0; j < 4; j ++)
					matLD[i][j] = (long double) genrand_real3();
			// Initialize the constants
			struct constants_t *consts_logmetrics = malloc(sizeof(struct constants_t));
			InitConstants(consts_logmetrics);

			// Initialize an array to hold the output metric values
			double *metvals = malloc(sizeof(double) * nmetrics);

			// Call function to compute the metric values.
			ComputeMetrics(metvals, nmetrics, metrics, matLD, chname, consts_logmetrics);

			// Print metric values
			for (i = 0; i < nmetrics; i ++)
				printf("%s = %g\n", metrics[i], metvals[i]);

			// Free memory
			free(chname);
			for (i = 0; i < nmetrics; i ++)
				free(metrics[i]);
			free(metrics);
			free(metvals);
			FreeConstants(consts_logmetrics);
			free(consts_logmetrics);
		}
	}

	if (strncmp(file, "memory", 6) == 0){
		// Testing functions in memory.c
		// Create a pointer to simulation structure.
		int nphys = 7, nenc = 1;
		struct simul_t *sim = malloc(sizeof(struct simul_t));
		if (strncmp(func, "AllocSimParams", 14) == 0){
			printf("Function: AllocSimParams.\n");
			AllocSimParams(sim, nphys, nenc);
		}
		if (strncmp(func, "AllocSimParamsQECC", 18) == 0){
			printf("Function: AllocSimParamsQECC.\n");
			AllocSimParamsQECC(sim, nphys, nenc);
		}
		if (strncmp(func, "FreeSimParams", 13) == 0){
			printf("Function: FreeSimParams.\n");
			FreeSimParams(sim, nphys, nenc);
		}
		if (strncmp(func, "FreeSimParamsQECC", 17) == 0){
			printf("Function: FreeSimParamsQECC.\n");
			FreeSimParamsQECC(sim, nphys, nenc);
		}
	}

	if (strncmp(file, "utils", 4) == 0){
		// Test the functions in the rand.c file.
		printf("Testing functions to perform simple arethematic operations.\n");
		if (strncmp(func, "Divide", 6) == 0){
			long double numerator, denominator, naive, ours, scale = 30;
			int t, trials = 10;
			for (t = 0; t < trials; t ++){
				numerator = powl(-1, (long double) RandomRangeInt(0, 2)) * ldexpl((long double) genrand_real3(), (-1) * RandomRangeInt(0, scale));
				denominator = powl(-1, (long double) RandomRangeInt(0, 2)) * ldexpl((long double) genrand_real3(), (-1) * RandomRangeInt(0, scale));
				printf("%d). A = %.15Le, B = %.15Le.\n", t + 1, numerator, denominator);
				naive = numerator/denominator;
				ours = Divide(numerator, denominator);
				printf("Naive A/B = %.15Le and our A/B = %.15Le. Difference: %.15Le.\n", naive, ours, fabsl(naive - ours));
			}
		}
		if (strncmp(func, "Order", 5) == 0){
			// Test the function to compute the order of magnitude of a number.
			const int base = 10, max_expo = 50;
			long double num;
			int t, trials = 10;
			for (t = 0; t < trials; t ++){
				num = powl(-1, (long double) RandomRangeInt(0, 2)) * ldexpl((long double) genrand_real3(), (-1) * RandomRangeInt(0, max_expo));
				printf("%d). The order of magnitude of %.5Le is %d.\n", t + 1, num, OrderOfMagnitude(num, base));
			}
		}
	}

	if (strncmp(file, "rand", 4) == 0){
		// Test the functions in the rand.c file.
		printf("Testing functions to generate random patterns.\n");
		if (strncmp(func, "ShuffleInt", 10) == 0){
			int size = 10, nshuffles = 10;
			int *arr = malloc(sizeof(int) * size);
			for (i = 0; i < size/2; i ++)
				arr[i] = 1;
			for (i = size/2; i < size; i ++)
				arr[i] = 0;
			PrintIntArray1D(arr, "Array before shuffle", size);
			ShuffleInt(arr, size, nshuffles);
			printf("After applying %d shuffles.\n", nshuffles);
			PrintIntArray1D(arr, "Array after shuffle", size);
			free(arr);
		}
	}

	if (strncmp(file, "linalg", 6) == 0){
		if (strncmp(func, "Diagonalize", 11) == 0){
			int dim = 4;
			for (i = 0; i < 4; i ++)
				for (j = 0; j < 4; j ++)
					mat[i][j] = genrand_real3() + genrand_real3() * I;
			double complex *eigvals = malloc(sizeof(double complex) * dim);
			double complex **eigvecs = malloc(sizeof(double complex *) * dim);
			for (i = 0; i < 4; i ++)
				eigvecs[i] = malloc(sizeof(double complex) * dim);

			PrintComplexArray2D(mat, "M", dim, dim);
			// DiagonalizeD(mat, dim, eigvals, 1, eigvecs);
			ZEigH(mat, dim, eigvals, 1, eigvecs);
			PrintSpectralDecomposition(eigvals, eigvecs, "Before reconstruction", 4);

			printf("Check with Python:\n");
			PrintPythonComplexArray2D(mat, "M", dim, dim);
			
			// Free memory
			free(eigvals);
			for (i = 0; i < 4; i ++)
				free(eigvecs[i]);
			free(eigvecs);
		}
		if (strncmp(func, "BinaryDot", 9) == 0){
			int max = 64;
			int a = RandomRangeInt(0, max);
			int *parity = malloc(max * sizeof(int));
			printf("Testing the binary dot product between %d and all number form 0 to %d.\n", a, max - 1);
			for (i = 0; i < max; i ++)
				parity[i] = BinaryDot(i, a);
			PrintIntArray1D(parity, "Parities", max);
			printf("Number of non-orthogonal numbers = %d\n", SumInt(parity, max));
			free(parity);
		}
		if (strncmp(func, "FixPositivity", 13) == 0){
			// Test the function that, when given M, computes M' such that M' has the positive eigenspectrum of M.
			double complex **nonpsd_mat = malloc(sizeof(double complex *) * 4);
			double complex **psd_mat = malloc(sizeof(double complex *) * 4);
			double complex **herm_mat = malloc(sizeof(double complex *) * 4);
			for (i = 0; i < 4; i ++){
				nonpsd_mat[i] = malloc(sizeof(double complex) * 4);
				psd_mat[i] = malloc(sizeof(double complex) * 4);
				herm_mat[i] = malloc(sizeof(double complex) * 4);
			}
			const int trials = 1, scale = 20;
			const double violation = 1E-3;
			double real_part, imag_part;
			int t, success;
			for (t = 0; t < trials; t ++){
				printf("%d).\n", t + 1);
				for (i = 0; i < 4; i ++){
					for (j = 0; j < 4; j ++){
						real_part = pow(-1, RandomRangeInt(0, 2)) * ldexp(genrand_real3(), (-1) * RandomRangeInt(scale - 4, scale));
						imag_part = pow(-1, RandomRangeInt(0, 2)) * ldexp(genrand_real3(), (-1) * RandomRangeInt(scale - 4, scale));
						nonpsd_mat[i][j] = real_part + I * imag_part;
						psd_mat[i][j] = 0;
						herm_mat[i][j] = 0;
					}
				}
				// Force the complex non-psd matrix to be Hermitian.
				// M -> M + M^\dag
				for (i = 0; i < 4; i ++)
					for (j = 0; j < 4; j ++)
						herm_mat[i][j] = nonpsd_mat[i][j] + conj(nonpsd_mat[j][i]);

				printf("Non-Positive matrix M\n");
				PrintComplexArray2D(herm_mat, "M", 4, 4);
				PrintPythonComplexArray2D(herm_mat, "M", 4, 4);
				success = FixPositivity(herm_mat, psd_mat, 4, violation);
				if (success == 1){
					printf("Positive semidefinite matrix M_+\n");
					PrintComplexArray2D(psd_mat, "M_+", 4, 4);
					PrintPythonComplexArray2D(psd_mat, "M_+", 4, 4);
					printf("||M - M_+||_2 = %.5e.\n", ZFroNorm(herm_mat, psd_mat, 4, 4));
				}
				else
					printf("Matrix is not even approximately CP, within a tolerance of %.2e.\n", violation);
				printf("xxxxxxxxxx\n");
			}

			// Free memory
			for (i = 0; i < 4; i ++){
				free(nonpsd_mat[i]);
				free(herm_mat[i]);
				free(psd_mat[i]);
			}
			free(nonpsd_mat);
			free(herm_mat);
			free(psd_mat);
		}
	}

	if (strncmp(file, "benchmark", 9) == 0){
		// Testing the entire benchmarking functionality.
		printf("Testing the entire benchmarking functionality.\n");
		// The array inputs for this function's test are in the folder: ./../input/debug_test/
		if (strncmp(func, "Benchmark", 9) == 0){
			int nlevels = 3;

			// ===
			int *nkd = malloc(nlevels * 3 * sizeof(int));
			LoadIntArray1D(nkd, "./../chflow/input/debug_testing/nkd.txt", nlevels * 3);
			PrintIntArray1D(nkd, "nkd", nlevels * 3);
			// ===

			// ===
			int n = nkd[0], k = nkd[1];
			int nstabs = (int) pow(2, (n - k)), nlogs = (int) pow(4, k);
			printf("nstabs = %d, nlogs = %d.\n", nstabs, nlogs);
			// ===

			// ===
			int *SS = malloc(nlevels * nstabs * nstabs * sizeof(int));
			LoadIntArray1D(SS, "./../chflow/input/debug_testing/SS.txt", nlevels * nstabs * nstabs);
			printf("_/ Loaded SS.\n");
			// PrintIntArray1D(SS, "SS", nlevels * nstabs * nstabs);
			// ===

			// ===
			int *normalizer = malloc(nlevels * nlogs * nstabs * n * sizeof(int));
			LoadIntArray1D(normalizer, "./../chflow/input/debug_testing/normalizer.txt", nlevels * nlogs * nstabs * n);
			printf("_/ Loaded normalizer.\n");
			// PrintIntArray1D(normalizer, "normalizer", nlevels * nlogs * nstabs * n);
			// ===

			// ===
			double *normphases_real = malloc(nlevels * nlogs * nstabs * sizeof(double));
			LoadDoubleArray1D(normphases_real, "./../chflow/input/debug_testing/normphases_real.txt", nlevels * nlogs * nstabs);
			printf("_/ Loaded normphases_real.\n");
			// PrintDoubleArray1D(normphases_real, "normphases_real", nlevels * nlogs * nstabs);
			// ===

			// ===
			double *normphases_imag = malloc(nlevels * nlogs * nstabs * sizeof(double));
			LoadDoubleArray1D(normphases_imag, "./../chflow/input/debug_testing/normphases_imag.txt", nlevels * nlogs * nstabs);
			printf("_/ Loaded normphases_imag.\n");
			// PrintDoubleArray1D(normphases_imag, "normphases_imag", nlevels * nlogs * nstabs);
			// ===

			// ===
			char *chname = malloc(100 * sizeof(char));
			sprintf(chname, "cptp");
			// ===

			// ===
			int iscorr = 3;
			// ===

			// ===
			int nparams = 0;
			if (iscorr == 0)
				nparams = nlogs * nlogs;
			else if (iscorr == 1)
				nparams = nstabs * nlogs;
			else if (iscorr == 2)
				nparams = n * nlogs * nlogs;
			else
				nparams = nstabs * nlogs * nstabs * nlogs;
			// ===

			// ===
			// printf("nparams = %d\n", nparams);
			double *physical = malloc(nparams * sizeof(double));
			LoadDoubleArray1D(physical, "./../chflow/input/debug_testing/physical.txt", nparams);
			printf("Loaded physical.\n");
			// PrintDoubleArray1D(physical, "physical", nparams);
			// ===

			// ===
			int nmetrics = 1;
			char **metrics = malloc(nmetrics * sizeof(char *));
			metrics[0] = malloc(100 * sizeof(char));
			sprintf(metrics[0], "infid");
			printf("metrics[0]: %s.\n", metrics[0]);
			// ===

			// ===
			int *decoders = malloc(nlevels * sizeof(int));
			for (i = 0; i < nlevels; i ++)
				decoders[i] = 1;
			PrintIntArray1D(decoders, "Decoders", nlevels);
			int *dclookups = malloc(nlevels * nstabs * sizeof(int));
			LoadIntArray1D(dclookups, "./../chflow/input/debug_testing/lookup.txt", nlevels * nstabs);
			int *operators_LST = malloc(nlevels * nstabs * nlogs * nstabs * n * sizeof(int));
			LoadIntArray1D(operators_LST, "./../chflow/input/debug_testing/lst.txt", nlevels * nstabs * nlogs * nstabs * n);
			double *mpinfo = malloc((int)pow(4, n) * sizeof(double));
			LoadDoubleArray1D(mpinfo, "./../chflow/input/debug_testing/mpinfo.txt", (int)pow(4, n));
			// PrintDoubleArray1D(mpinfo, "Message passing initialize", (int)pow(4, n));
			// ===

			// ===
			int hybrid = 0;
			int *decoderbins = NULL;
			int *ndecoderbins = NULL;
			// ===

			// ===
			int frame = 0;
			// ===

			// ===
			int nbreaks = 1;
			// ===

			// ===
			long *stats = malloc(nbreaks * sizeof(long));
			stats[0] = 10000;
			printf("Stats = %ld\n", stats[0]);
			// ===

			// ===
			int nbins = 1;
			int maxbin = 1;
			// ===

			// ===
			int importance = 1;
			printf("Importance = %d\n", importance);
			// ===

			// ===
			double *refchan = NULL;
			// ===

			// ===
			int rc = 0;
			// ===

			//
			double infidelity = 0.020499323636302913;
			//

			// Calling the Benchmark function
			int trials = 1;
			for (i = 0; i < trials; i ++){
				printf("Calling the Benchmark function for the %d time with %d levels and rc = %d.\n", i+1, nlevels, rc);
				Benchmark(nlevels, nkd, SS, normalizer, normphases_real, normphases_imag, chname, iscorr, physical, rc, nmetrics, metrics, decoders, dclookups, mpinfo, operators_LST, hybrid, decoderbins, ndecoderbins, frame, nbreaks, stats, nbins, maxbin, importance, refchan, infidelity);
			}
			// Free memory
			free(nkd);
			free(SS);
			free(normalizer);
			free(normphases_real);
			free(normphases_imag);
			free(chname);
			free(physical);
			free(metrics[0]);
			free(metrics);
			free(stats);
			free(decoders);
			free(dclookups);
			free(mpinfo);
			free(operators_LST);
		}
	}

	// Free memory
	for (i = 0; i < 4; i ++)
		free(mat[i]);
	free(mat);
	for (i = 0; i < 4; i ++)
		free(matLD[i]);
	free(matLD);
	free(file);
	free(func);

	// Output the total run time of the test.
	clock_t end = clock();
	double runtime = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("***********\n");
	printf("All testing done in %d seconds.\n", (int) runtime);
	printf("***********\n");

	return 0;
}
