#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include "reps.h"

void ChoiToProcess(double **process, double complex **choi, double complex ***pauli)
{
	// Convert from the Choi matrix to the process matrix, of a quantum channel.
	// Gamma[a,b] = Trace( Choi * (Pb \otimes Pa^T) ).
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

void ProcessToChoi(double **process, double complex **choi, int nlogs, double complex ***pauli){
	// Convert from the process matrix to the Choi matrix.
	// J = 1/K * sum_P (E(P) o P^T)
	// where K is the dimension of the Hilbert space, P runs over Pauli operators
	// and E is the error channel. E(P) = 1/K * sum_Q G_(P,Q) Q where G(P, Q) =
	// Tr(E(P).Q) is the process matrix and Q runs over Pauli operators. hence we
	// have: J = 1/K sum_(P, Q) G_(P,Q) Q o P^T
	// printf("Function: ProcessToChoi.\n");
	int v, i, j, p, q;
	for (v = 0; v < nlogs * nlogs; v++)
		choi[v / nlogs][v % nlogs] = 0;
	for (v = 0; v < nlogs * nlogs * nlogs * nlogs; v++) {
		p = v % nlogs;
		q = (v / nlogs) % nlogs;
		j = (v / (nlogs * nlogs)) % nlogs;
		i = (v / (nlogs * nlogs * nlogs)) % nlogs;
		choi[i][j] += (double complex)(0.25 * process[p][q]) * (double complex) (pauli[q][i / 2][j / 2]) * (double complex) (pauli[p][j % 2][i % 2]);
	}
	// PrintComplexArray2D(choi, "Choi", nlogs, nlogs);
}
