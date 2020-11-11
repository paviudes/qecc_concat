#ifndef SAMPLING_H
#define SAMPLING_H

/*
	Computing an exponent of the true distribution such that outlier events whose probability is of the order of the average infidelity, is increased to a threshold: \lambda..
	1. Compute an upperbound to 1 - F, for the given level and physical infidelity.
	2. The exponent is give by: k = log(\lambda)/log(1 - F).
*/
extern double SetExponent(double phy_infid, int dist, int nphys, int level);

/*
	Setting outlier syndrome probabilities.
	Upper and lower limits for the probability of the outlier syndromes.
	We will assume that the probability of the outlier syndromes is not more than p^d/2 and not less than 80% of p^d/2
	For a physical noise process whose Pauli transfer matrix is G, we will define p = 0.5 + 0.5 * (4 - tr(G))/4.
	Additionally, we want to make sure that 0.5 <= p <= 1. This is safe for the importance sampler since p ~ 0 will lead to an indefinite search in PowerSearch(...) in sampling.c.
	We will follow the definition of infidelity in eq. 5.16 of https://arxiv.org/abs/1109.6887.pdf.		
*/
extern void SetOutlierProbs(double phy_infid, int dist, int level, double *outlierprobs);

/*
	Given an exponent k and a probability distribution P(s), construct a new probability distribution Q(s) where
	Q(s) = P(s)^k/(sum_s P(s)^k), i.e, a normalized power-law scaled version of P(s).
	If k = 0 (i.e, less than 10E-5) then simply set the probability distribution to be flat.
*/
extern void ConstructImportanceDistribution(double* truedist, double *impdist, int nelems, double expo);

/*
	Given a probability distribution, construct its cumulative distribution.
*/
extern void ConstructCumulative(double* dist, double *cumul, int nelems);

/*
	Sample a discrete probability distribution given its cumulative distribution.
	Draw a uniform random number, u, in [0,1]. Determine the interval of the cumulative distribution in which u lies.
	http://ieeexplore.ieee.org/document/92917/
	Additions: if frozen = -1, continue with the sampling as described above.
*/
extern int SampleCumulative(double *cumulative, int size);

/*
	Determine is a number is within, below or above a window.
	If it is within, return 0, if it is above, return 1 and if below, return -1.
*/
extern int WhereInWindow(double number, double *window);

/*
	Search for an exponent k such that according to normalized distribution of P(s)^k, the probability of isincluded errors is within a desired window.
	the array searchin provides the exponent k in that: k = (searchin[0] + searchin[1])/2.
	One of the three cases can occur for the distribution P(s)^k where k = (searchin[0] + searchin[1])/2.
		1. Probability of isincluded errors is below the window -- in this case, we return the value of the function on a new searchin wondow, given by: [searchin[0], k]. [Recursion]
		2. Probability of isincluded errors is within the window -- in this case we stop after returning k.
		3. Probability of isincluded errors is above the window -- in this case, we return the value of the function on a new searchin wondow, given by: [k, searchin[1]]. [Recursion]
	The weight-1 errors are X_i, Z_i, Y_i for i = 1 to 7 and their syndromes are: 56, 24, 40, 8, 48, 16, 32, 7, 3, 5, 1, 6, 2, 4, 63, 54, 45, 36, 27, 18, 9, respectively.
*/
extern double PowerSearch(double *dist, int size, double window[2], double searchin[2]);

#endif /* SAMPLING_H */