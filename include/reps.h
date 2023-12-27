#ifndef REPS_H
#define REPS_H

#include <complex.h>
#include "constants.h"

/*
	Convert from the Choi matrix to the process matrix, of a quantum channel.
	CHI[a,b] = Trace( Choi * (Pb \otimes Pa^T) ).
*/
extern void ChoiToProcess(double **process, double complex **choi, int nlogs, double complex ***pauli);

/*
	Convert from the process matrix to the Choi matrix.
	J = 1/K * sum_P (E(P) o P^T)
	where K is the dimension of the Hilbert space, P runs over Pauli operators
	and E is the error channel. E(P) = 1/K * sum_Q G_(P,Q) Q where G(P, Q) =
	Tr(E(P).Q) is the process matrix and Q runs over Pauli operators. hence we
	have: J = 1/K sum_(P, Q) G_(P,Q) Q o P^T
*/
extern void ProcessToChoi(double **process, double complex **choi, int nlogs, double complex ***pauli);

#endif /* REPS_H */
