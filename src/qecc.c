#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "rand.h"
#include "checks.h" // only for testing purposes.
#include "constants.h"
#include "linalg.h"
#include "memory.h"
#include "printfuns.h"
#include "qecc.h"

void InitQECC(struct qecc_t *qecc)
{
	// Initialize the elements of a quantum error correcting code.
	qecc->nlogs = (int)pow(4, (double)(qecc->K));
	qecc->nstabs = (int)pow(2, (double)(qecc->N - qecc->K));

	// printf("Function: InitQECC: nlogs = %d, nstabs = %d.\n", qecc->nlogs,
	// qecc->nstabs);

	int s;
	qecc->projector = malloc(qecc->nstabs * sizeof(int *));
	for (s = 0; s < qecc->nstabs; s++)
		qecc->projector[s] = malloc(qecc->nstabs * sizeof(int));

	int i;
	qecc->action = malloc(qecc->nlogs * sizeof(int **));
	for (i = 0; i < qecc->nlogs; i++)
	{
		(qecc->action)[i] = malloc(qecc->nstabs * sizeof(int *));
		for (s = 0; s < qecc->nstabs; s++)
			(qecc->action)[i][s] = malloc(qecc->N * sizeof(int));
	}

	qecc->phases = malloc(qecc->nlogs * sizeof(double complex *));
	for (i = 0; i < qecc->nlogs; i++)
		(qecc->phases)[i] = malloc(qecc->nstabs * sizeof(double complex));

	qecc->dclookup = malloc(qecc->nstabs * sizeof(int));
	for (s = 0; s < qecc->nstabs; s++)
		qecc->dclookup[s] = 0;

	qecc->dcknowledge = malloc(qecc->nstabs * qecc->nlogs * sizeof(double));
	for (s = 0; s < qecc->nstabs * qecc->nlogs; s++)
		qecc->dcknowledge[s] = 0.0;

	// printf("Done InitQECC.\n");
}

void FreeQECC(struct qecc_t *qecc)
{
	// Free the memory assigned to a quantum error correcting code.
	int s;
	for (s = 0; s < qecc->nstabs; s++)
		free((qecc->projector)[s]);
	free(qecc->projector);

	int i;
	for (i = 0; i < qecc->nlogs; i++)
	{
		for (s = 0; s < qecc->nstabs; s++)
			free((qecc->action)[i][s]);
		free((qecc->action)[i]);
	}
	free(qecc->action);

	for (i = 0; i < qecc->nlogs; i++)
		free((qecc->phases)[i]);
	free(qecc->phases);

	free(qecc->dclookup);
	free(qecc->dcknowledge);
}

void ChoiToProcess(double **process, double complex **choi, double complex ***pauli)
{
	// Convert from the Choi matrix to the process matrix, of a quantum channel.
	// CHI[a,b] = Trace( Choi * (Pb \otimes Pa^T) ).
	int v, i, j, k;
	double complex contribution;
	// #pragma ivdep
	for (v = 0; v < 16; v++)
	{
		j = v % 4;
		i = v / 4;
		contribution = 0;
		for (k = 0; k < 16; k++)
			contribution += choi[k / 4][k % 4] * pauli[j][(k % 4) / 2][(k / 4) / 2] * pauli[i][(k / 4) % 2][(k % 4) % 2];
		process[i][j] = creal(contribution);
	}
}

void GetFullProcessMatrix(struct qecc_t *qecc, struct simul_t *sim, int isPauli)
{
	// For each pair of logical operators, we only need the entries of the Chi
	// matrix that correspond to Pauli operators from different logical classes.
	// We construct the sections of the Chi matrices that correspond to rows and
	// columns in the logical classes. printf("Function: GetFullProcessMatrix for
	// isPauli = %d.\n", isPauli);
	int i, j, k, l, q;
	double contribution = 1;
	if (isPauli == 0)
	{
		for (i = 0; i < qecc->nlogs; i++)
		{
			for (j = 0; j < qecc->nlogs; j++)
			{
				for (k = 0; k < qecc->nstabs; k++)
				{
					for (l = 0; l < qecc->nstabs; l++)
					{
						contribution = 1;
						for (q = 0; q < qecc->N; q++)
							contribution *= (sim->virtchan)[q][(qecc->action)[i][k][q]][(qecc->action)[j][l][q]];
						(sim->process)[i][j][k][l] = creal((qecc->phases)[i][k] * (qecc->phases)[j][l]) * contribution;
					}
				}
			}
		}
	}
	else
	{
		// For a Pauli channel, the process matrix is diagonal.
		for (i = 0; i < qecc->nlogs; i++)
		{
			for (j = 0; j < qecc->nstabs; j++)
			{
				contribution = 1;
				for (q = 0; q < qecc->N; q++)
					contribution *= (sim->virtchan)[q][(qecc->action)[i][j][q]][(qecc->action)[i][j][q]];
				(sim->process)[i][i][j][j] = creal((qecc->phases)[i][j] * (qecc->phases)[i][j]) * contribution;
			}
		}
		// PrintDoubleArrayDiag((sim->process)[1][1], "process[1][1]", qecc->nstabs);
	}
	if (sim->rc == 1)
	{
		// PrintIntArray1D(sim->rcpauli, "rcpauli", qecc->nstabs);
		for (i = 0; i < qecc->nlogs; i++)
		{
			for (j = 0; j < qecc->nlogs; j++)
			{
				for (k = 0; k < qecc->nstabs; k++)
				{
					for (l = 0; l < qecc->nstabs; l++)
					{
						(sim->process)[i][j][k][l] *= pow(-1, (sim->rcpauli)[k]);
					}
				}
			}
		}
		// PrintDoubleArray2D((sim->process)[0][0], "process[0][0]", qecc->nstabs, qecc->nstabs);
	}
	// printf("Done computing full process matrix: process[0][0][0][0] = %g.\n", (sim->process)[0][0][0][0]);
}


void ComputeSyndromeProbability(int synd, struct qecc_t *qecc, struct simul_t *sim, int isPauli)
{
	// Compute the probability of all the syndromes in the qecc code, for the
	// given error channel and input state. Sample a syndrome from the resulting
	// probability distribution of the syndromes. Probability of a syndrome s,
	// denoted by P(s) is given by the following expression. P(s) = 1/2^(n-k) *
	// sum_(i,j: P_i and P_j are stabilizers) CHI[i,j] * (-1)^sign(P_j).
	// Initialize syndrome probabilities
	if (isPauli == 0)
		(sim->syndprobs)[synd] = SumDotInt((sim->process)[0][0], (qecc->projector)[synd], qecc->nstabs, qecc->nstabs, qecc->nstabs) / (double)(qecc->nstabs);
	else
		(sim->syndprobs)[synd] = DiagGDotIntV((sim->process)[0][0], (qecc->projector)[synd], qecc->nstabs, qecc->nstabs, qecc->nstabs) / (double)(qecc->nstabs);
}


void ComputeSyndromeDistribution(struct qecc_t *qecc, struct simul_t *sim, int isPauli)
{
	// Compute the probability of all the syndromes in the qecc code, for the
	// given error channel and input state. Sample a syndrome from the resulting
	// probability distribution of the syndromes. Probability of a syndrome s,
	// denoted by P(s) is given by the following expression. P(s) = 1/2^(n-k) *
	// sum_(i,j: P_i and P_j are stabilizers) CHI[i,j] * (-1)^sign(P_j).
	int s;
	// Initialize syndrome probabilities
	for (s = 0; s < qecc->nstabs; s++)
		ComputeSyndromeProbability(s, qecc, sim, isPauli);
	// Construct the cumulative distribution
	(sim->cumulative)[0] = (sim->syndprobs)[0];
	for (s = 1; s < qecc->nstabs; s++)
		(sim->cumulative)[s] = (sim->cumulative)[s - 1] + (sim->syndprobs)[s];
	// PrintDoubleArray1D((sim->syndprobs), "Syndrome distribution", qecc->nstabs);
	// PrintDoubleArray1D((sim->cumulative), "Cumulative Syndrome distribution", qecc->nstabs);
}


/*
void MLDecoder(struct qecc_t *qecc, struct simul_t *sim, struct constants_t
*consts, int currentframe, int isPauli){
	   // Perform maximum likelihood decoding.
	   // Compute the probabilities of the logical classes, considitioned on a
particular syndrome.
	   // The ML Decoder picks the logical error which belongs to the class
that has the maximum probability.
	   // The probability of a logical class is P(L|s) = Tr( L r L . Ts PI_s
E(r) PI_s Ts )/P(s) which can be simplified to
	   // P(L|s) = 1/P(s) * sum_(u: Paulis) sum_(i: P_i is in the [u] logical
class) sum_(j: Pj is in the [L u L] logical class) CHI[i,j] * (-1)^(P_j).
	   // inputs: nqecc, kqecc, chi, algebra (conjugations).
	   // printf("Function: MLDecoder\n");
	   const double atol = 10E-10;
	   int v, i, j, u, l, s;
	   double prob, maxprob;
	   for (s = 0; s < qecc->nstabs; s ++){
			 if ((sim->syndprobs)[s] > atol){
				    (sim->corrections)[s] = 0;
				    maxprob = 0;
				    for (l = 0; l < currentframe; l ++){
						  prob = 0;
						  if (isPauli == 0)
								for (u = 0; u < qecc->nlogs; u ++)
									   prob +=
(consts->algebra)[1][l][u] *
SumDotInt((sim->process)[u][(consts->algebra)[0][l][u]], (qecc->projector)[s],
qecc->nstabs, qecc->nstabs, qecc->nstabs); else for (u = 0; u < qecc->nlogs; u
++) if ((consts->algebra)[0][l][u] == u) prob += (consts->algebra)[1][l][u] *
DiagGDotIntV((sim->process)[u][u], (qecc->projector)[s], qecc->nstabs,
qecc->nstabs, qecc->nstabs); printf("s = %d, l = %d, prob = %g\n", s, l,
prob/(qecc->nlogs * qecc->nstabs * (sim->syndprobs)[s])); if (prob > maxprob){
								(sim->corrections)[s] = l;
								maxprob = prob;
						  }
				    }
				    printf("\n");
			 }
	   }
	   PrintIntArray1D(sim->corrections, "sim->corrections", qecc->nstabs);
}
*/


void MLDecodeSyndrome(int synd, struct qecc_t *qecc, struct simul_t *sim, struct constants_t *consts, int currentframe, int isPauli, int dcalg)
{
	// Perform maximum likelihood decoding.
	// Compute the probabilities of the logical classes, considitioned on a
	// particular syndrome. The ML Decoder picks the logical error which belongs
	// to the class that has the maximum probability. The probability of a logical
	// class is P(L|s) = Tr( L r L . Ts PI_s E(r) PI_s Ts )/P(s) which can be
	// simplified to P(L|s) = 1/P(s) * sum_(u: Paulis) sum_(i: P_i is in the [u]
	// logical class) sum_(j: Pj is in the [L u L] logical class) CHI[i,j] *
	// (-1)^(P_j). inputs: nqecc, kqecc, chi, algebra (conjugations).
	// printf("Function: MLDecodeSyndrome %d, dcalg = %d, currentframe = %d\n", synd, dcalg, currentframe);
	const double atol = 10E-16;
	int i, j, u, l;
	double prob, maxprob, contrib;
	(sim->corrections)[synd] = 0;
	if ((sim->syndprobs)[synd] > atol)
	{
		maxprob = 0;
		for (l = 0; l < currentframe; l++)
		{
			prob = 0;
			if (isPauli == 0)
			{
				for (u = 0; u < qecc->nlogs; u++)
				{
					contrib = 0;
					for (i = 0; i < qecc->nstabs; i++)
						for (j = 0; j < qecc->nstabs; j++)
							contrib += (sim->process)[u][(consts->algebra)[0][l][u]][i][j] * (qecc->projector)[synd][j];
					prob += (consts->algebra)[1][l][u] * contrib;
				}
			}
			else
			{
				if (dcalg == 0){
					for (u = 0; u < qecc->nlogs; u++)
					{
						contrib = 0;
						for (i = 0; i < qecc->nstabs; i++)
							contrib = contrib + (sim->process)[u][u][i][i] * (qecc->projector)[synd][i];
						prob = prob + (consts->algebra)[1][l][u] * contrib;
					}
				}
				else{
					// Case of the partial knowledge decoder, i.e., dcalg == 2.
					for (u = 0; u < qecc->nlogs; u++)
					{
						contrib = 0;
						for (i = 0; i < qecc->nstabs; i++)
							contrib = contrib + (qecc->dcknowledge)[u * qecc->nstabs + i] * (qecc->projector)[synd][i];
						prob = prob + (consts->algebra)[1][l][u] * contrib;
					}
				}
			}
			if (prob > maxprob)
			{
				(sim->corrections)[synd] = l;
				maxprob = prob;
			}
		}
	}
}


void MLDecoder(struct qecc_t *qecc, struct simul_t *sim, struct constants_t *consts, int dcalg, int currentframe, int isPauli)
{
	// Perform maximum likelihood decoding.
	// Compute the probabilities of the logical classes, considitioned on a
	// particular syndrome. The ML Decoder picks the logical error which belongs
	// to the class that has the maximum probability. The probability of a logical
	// class is P(L|s) = Tr( L r L . Ts PI_s E(r) PI_s Ts )/P(s) which can be
	// simplified to P(L|s) = 1/P(s) * sum_(u: Paulis) sum_(i: P_i is in the [u]
	// logical class) sum_(j: Pj is in the [L u L] logical class) CHI[i,j] *
	// (-1)^(P_j). inputs: nqecc, kqecc, chi, algebra (conjugations).
	// printf("Function: MLDecoder\n");
	int s;
	for (s = 0; s < qecc->nstabs; s++)
	{
		if (dcalg == 0)
			MLDecodeSyndrome(s, qecc, sim, consts, currentframe, isPauli, dcalg);
		// else if (dcalg == 2) // the partial knowledge decoder only accesses the diagonal entries.
		// 	MLDecodeSyndrome(s, qecc, sim, consts, currentframe, 1, dcalg);
		else
			(sim->corrections)[s] = (qecc->dclookup)[s];
	}
	// PrintIntArray1D(sim->corrections, "Corrections after decoding", qecc->nstabs);
	// printf("dcalg = %d.\n", dcalg);
}


void EffChanSynd(int synd, struct qecc_t *qecc, struct simul_t *sim, struct constants_t *consts, int isPauli){
	// Compute the effective channel given a syndrome.
	// The effective channel in the Pauli Liouville representation is given by:
	// G_{L,L'} = \sum_S,S' G_{LS, L* L' L* S'}/P(s)
	// where L* is the correction applied by the decoder.
	// printf("Function: EffChanSynd\n");
	int l, lp, s, sp;
	int f1, f2, f3;
	// Initialization
	for (l = 0; l < qecc->nlogs; l ++){
		for (lp = 0; lp < qecc->nlogs; lp ++){
			(sim->effprocess)[synd][l][lp] = 0;
		}
	}
	// printf("Initialization done.\n");
	if (isPauli == 0){
		for (l = 0; l < qecc->nlogs; l ++){
			for (lp = 0; lp < qecc->nlogs; lp ++){
				for (s = 0; s < qecc->nstabs; s ++){
					for (sp = 0; sp < qecc->nstabs; sp ++){
						f1 = (qecc->projector)[synd][sp];
						f2 = (consts->algebra)[1][(sim->corrections)[synd]][lp];
						f3 = (consts->algebra)[0][(sim->corrections)[synd]][lp];
						(sim->effprocess)[synd][l][lp] += f1 * f2 * (sim->process)[l][f3][s][sp];
					}
				}
			}
		}
	}
	else{
		for (l = 0; l < qecc->nlogs; l ++){
			for (s = 0; s < qecc->nstabs; s ++){
				// printf("l = %d, s = %d\n",l,s);
				f1 = (qecc->projector)[synd][s];
				f2 = (consts->algebra)[1][(sim->corrections)[synd]][l];
				(sim->effprocess)[synd][l][l] += f1 * f2 * (sim->process)[l][l][s][s];
			}
		}
	}
	// printf("Population done.\n");
	// Normalization
	for (l = 0; l < qecc->nlogs; l ++){
		for (lp = 0; lp < qecc->nlogs; lp ++){
			(sim->effprocess)[synd][l][lp] /= (pow(2, qecc->N - qecc->K) * (sim->syndprobs)[synd]);
		}
	}
	// Consistency checks.
	// 1. Check G[0,0] = 1
	// printf("s = %d, P(s) = %.10f\nG[0][0] = %.18f\n", synd, (sim->syndprobs)[synd], (sim->effprocess)[synd][0][0]);
	// // 2. First column of G should be all zeros
	// for (l = 1; l < qecc->nlogs; l ++){
	// 	printf("G[%d][0] = %.18f\n", l, (sim->effprocess)[synd][l][0]);
	// }
	// printf("-------\n");
}


void ComputeEffectiveChannels(struct qecc_t *qecc, struct simul_t *sim, struct constants_t *consts, int isPauli){
	/*
		Compute effective channels for all syndromes.
		degeneracies = {l: [all syndromes for which l is the correction]}
		for l in degeneracies:
			// Compute the effective channel for degeneracies[l][0],
			// and store it as E.
			for i in range(1, len(degeneracies[l][0])):
				// let s = degeneracies[l][i]
				// effective channel for the syndrome s is E.
	*/
	int s;
	for (s = 0; s < qecc->nstabs; s++)
		EffChanSynd(s, qecc, sim, consts, isPauli);
}



void SetFullProcessMatrix(struct qecc_t *qecc, struct simul_t *sim, double *process, int isPauli)
{
	// In cases where the full process matrix needs to be explicity specified,
	// load the structure attribute with the specified matrix. In the case of
	// Pauli channels, the process matrix input must be a 1 x 4^n, each column
	// specifying a unique Pauli error probability, where the errors are ordered
	// as SLT.
	int l1, l2, s1, s2;
	int nlogs = qecc->nlogs, nstabs = qecc->nstabs;
	// PrintDoubleArray1D((sim->physical), "sim->physical", nlogs * nlogs * nstabs * nstabs);
	if (isPauli == 0)
		for (l1 = 0; l1 < nlogs; l1 ++)
			for (l2 = 0; l2 < nlogs; l2 ++)
				for (s1 = 0; s1 < nstabs; s1 ++)
					for (s2 = 0; s2 < nstabs; s2 ++)
						(sim->process)[l1][l2][s1][s2] = process[l1 * nlogs * nstabs * nstabs + s1 * nlogs * nstabs + l2 * nstabs + s2];
	else
		for (l1 = 0; l1 < nlogs; l1 ++)
			for (s1 = 0; s1 < nstabs; s1 ++)
				(sim->process)[l1][l1][s1][s1] = process[l1 * nstabs + s1];
	// printf("Full process matrix set for isPauli = %d.\n", isPauli);
}

void SingleShotErrorCorrection(int isPauli, int iscorr, int dcalg, int frame, struct qecc_t *qecc, struct simul_t *sim, struct constants_t *consts)
{
	// Compute the effective logical channel, when error correction is applied over a set of input physical channels.
	// printf("Constructing the full process matrix\n");
	if ((iscorr == 0) || (iscorr == 2))
		GetFullProcessMatrix(qecc, sim, isPauli);
	// Compute the probabilities of all the syndromes.
	ComputeSyndromeDistribution(qecc, sim, isPauli);
	// Maximum Likelihood Decoding (MLD) -- For every syndrome, compute the
	// probabilities of the logical classes and pick the one that is most likely.
	// printf("Maximum likelihood decoding\n");
	MLDecoder(qecc, sim, consts, dcalg, frame, isPauli);
	// For every syndrome, apply the correction and compute the new effective channel.
	// printf("Computing effective channel\n");
	ComputeEffectiveChannels(qecc, sim, consts, isPauli);
}
