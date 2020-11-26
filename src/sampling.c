#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mt19937/mt19937ar.h"
// #include "printfuns.h" // only for testing
#include "utils.h"
#include "sampling.h"

void SetOutlierProbs(double phy_infid, int dist, int level, double *outlierprobs){
	// Setting outlier syndrome probabilities.
	// Upper and lower limits for the probability of the outlier syndromes.
	// We will assume that the probability of the outlier syndromes is not more than p^d/2 and not less than 80% of p^d/2
	// For a physical noise process whose Pauli transfer matrix is G, we will define p = 0.5 + 0.5 * (4 - tr(G))/4.
	// Additionally, we want to make sure that 0.5 <= p <= 1. This is safe for the importance sampler since p ~ 0 will lead to an indefinite search in PowerSearch(...) in sampling.c.
	// We will follow the definition of infidelity in eq. 5.16 of https://arxiv.org/abs/1109.6887.pdf.		
	double infidelity = pow(phy_infid, floor((pow(dist, level - 1) + 1)/2));
	// double infidelity = phy_infid;
	outlierprobs[1] = 0.99 * infidelity; // We are assuming that infidelity = alpha, here.
	// outlierprobs[1] = 0.2 * (1 - Min(3 * infidelity, 1 - 1E-5));
	outlierprobs[0] = 0.80 * outlierprobs[1];
}

double UpperBoundInfidelity(double p, int d, int n, int level){
	// Upper bound the infidelity by simply taking 1 - probabability of correctable errors.
	double sum_corr;
	int i, l;
	for (l = 0; l < level; ++l){
		sum_corr = 0;
		for(i = 0; i <= (d-1)/2; ++i)
			sum_corr += (double) Comb(n, i) * pow(p, i) * pow(1 - p, n - i);
		p = 1 - sum_corr;
	}
	return 1 - sum_corr;
}

double SetExponent(double phy_infid, int dist, int nphys, int level){
	// Computing an exponent of the true distribution such that outlier events whose probability is of the order of the average infidelity, is increased to a threshold: \lambda..
	// 1. Compute an upperbound to 1 - F, for the given level and physical infidelity.
	// 2. The exponent is give by: k = log(\lambda)/log(1 - F).
	double p = 1 - pow(1 - phy_infid, 1/(double) nphys);
	const double threshold = 0.01;
	double infidelity = UpperBoundInfidelity(p, dist, nphys, level);
	double expo = log(threshold)/log(infidelity);
	printf("L = %d, 1 - F = %.8f, exponent = %.3f\n", level, infidelity, expo);
	return expo;
}

void ConstructImportanceDistribution(double *truedist, double *impdist, int nelems, double expo){
	// Given an exponent k and a probability distribution P(s), construct a new probability distribution Q(s) where
	// Q(s) = P(s)^k/(sum_s P(s)^k), i.e, a normalized power-law scaled version of P(s).
	// If k = 0 (i.e, less than 10E-5) then simply set the probability distribution to be flat.
	int s;
	double atol = 10E-8;
	if (expo > atol)
		for (s = 0; s < nelems; s ++)
			impdist[s] = pow(truedist[s], expo);
	else
		for (s = 0; s < nelems; s ++)
			impdist[s] = 1;
	// Normalization
	double norm = 0.0;
	for (s = 0; s < nelems; s ++)
		norm = norm + impdist[s];
	// printf("Normalization = %g\n", norm);
	for (s = 0; s < nelems; s ++)
		impdist[s] = impdist[s]/norm;
	// PrintDoubleArray1D(impdist, "impdist", nelems);
}

void ConstructCumulative(double* dist, double *cumul, int nelems){
	// Given a probability distribution, construct its cumulative distribution.
	int s;
	cumul[0] = dist[0];
	for (s = 1; s < nelems; s ++)
		cumul[s] = cumul[s - 1] + dist[s];
}

int SampleCumulative(double *cumulative, int size){
	// Sample a discrete probability distribution given its cumulative distribution.
	// Draw a uniform random number, u, in [0,1]. Determine the interval of the cumulative distribution in which u lies.
	// http://ieeexplore.ieee.org/document/92917/
	// Additions: if frozen = -1, continue with the sampling as described above.
	int i = 0;
	double urand = genrand_real3();
	for (i = 0; i < size; i ++)
		if (urand < cumulative[i])
			return i;
	return (size - 1);
}

int WhereInWindow(double number, double window[2]){
	// Determine is a number is within, below or above a window.
	// If it is within, return 0, if it is above, return 1 and if below, return -1.
	if (number < window[0])
		return -1;
	else
		if (number > window[1])
			return 1;
	return 0;
}

double PowerSearch(double *dist, int size, double *window, double *searchin){
	// Search for an exponent k such that according to normalized distribution of P(s)^k, the probability of isincluded errors is within a desired window.
	// the array searchin provides the exponent k in that: k = (searchin[0] + searchin[1])/2.
	// One of the three cases can occur for the distribution P(s)^k where k = (searchin[0] + searchin[1])/2.
		// 1. Probability of isincluded errors is below the window -- in this case, we return the value of the function on a new searchin wondow, given by: [searchin[0], k]. [Recursion]
		// 2. Probability of isincluded errors is within the window -- in this case we stop after returning k.
		// 3. Probability of isincluded errors is above the window -- in this case, we return the value of the function on a new searchin wondow, given by: [k, searchin[1]]. [Recursion]
	// The weight-1 errors are X_i, Z_i, Y_i for i = 1 to 7 and their syndromes are: 56, 24, 40, 8, 48, 16, 32, 7, 3, 5, 1, 6, 2, 4, 63, 54, 45, 36, 27, 18, 9, respectively.
	// printf("Function: PowerSearch with distribution of size %d to bring P(0) in the window [%g, %g], with an exponent in [%g, %g].\n", size, window[0], window[1], searchin[0], searchin[1]);
	// PrintDoubleArray1D(dist, "P(s)", size);
	int i, position = 0, exit = 0;
	double pinc = 0, atol = 10E-5, exponent = (searchin[0] + searchin[1])/((double)2);
	double *powerdist = malloc(sizeof(double) * size);
	// incl is an indicator array for un-correctable errors.
	int *incl = malloc(sizeof(int) * size);
	incl[0] = 0;
	for (i = 1; i < size; i ++)
		incl[i] = 1;

	// If window is [1, 1], we will simply return 0 -- the importance distribution is flat.
	if ((window[0] == 1) && (window[1] == 1)){
		exponent = 0;
		exit = 1;
	}

	// If searchin is [0, 1], we will assume it is the first call.
	if ((searchin[0] == 0) && (searchin[1] == 1)){
		// If the probability of included errors is above than the window to begin with, then there is no need of importance sampling.
		pinc = 0;
		for (i = 0; i < size; i ++)
			if (incl[i] == 1)
				pinc += dist[i];
		// printf("pinc = %g.\n", pinc);
		if (pinc >= window[0]){
			exponent = 1;
			exit = 1;
		}
	}

	// printf("Function: PowerSearch for an exponent k in [%g, %g] such that %g <= 1 - P_k(s=0).\n", searchin[0], searchin[1], window[0]);

	while (exit == 0){
		if (fabs(searchin[0] - searchin[1]) < atol)
			exit = 1;
		else{
			ConstructImportanceDistribution(dist, powerdist, size, exponent);
			// PrintDoubleArray1D(powerdist, "Q(s)", size);
			for (i = 0; i < size; i ++)
				if (incl[i] == 1)
					pinc += powerdist[i];
			// printf("pinc = %g.\n", pinc);

			position = WhereInWindow(pinc, window);
			if (position == 0){
				// Free the local variables before returning
				free(incl);
				free(powerdist);
				return exponent;
			}

			else if (position == -1){
				// If pinc is too small, reduce the exponent, i.e, take a higher root.
				searchin[1] = exponent;
				exponent = (searchin[0] + searchin[1])/((double)2);
				// Free the local variables before returning
				free(incl);
				free(powerdist);
				return PowerSearch(dist, size, window, searchin);
			}
			else if (position == 1){
				// If pinc is too large, increase the exponent, (make the exponent closer to one) so that it goes closer to its real value.
				searchin[0] = exponent;
				exponent = (searchin[0] + searchin[1])/((double)2);
				// Free the local variables before returning
				free(incl);
				free(powerdist);
				return PowerSearch(dist, size, window, searchin);
			}
			else
				continue;
		}
	}

	// Free local variables before ending.
	free(incl);
	free(powerdist);

	return exponent;
}
