#ifndef LOGMETRICS_H
#define LOGMETRICS_H

#include <complex.h>
#include "constants.h"

/*
	Convert from the process matrix to the Choi matrix.
	J = 1/K * sum_P (E(P) o P^T)
	where K is the dimension of the Hilbert space, P runs over Pauli operators
	and E is the error channel. E(P) = 1/K * sum_Q G_(P,Q) Q where G(P, Q) =
	Tr(E(P).Q) is the process matrix and Q runs over Pauli operators. hence we
	have: J = 1/K sum_(P, Q) G_(P,Q) Q o P^T
*/
extern void ProcessToChoi(double **process, double complex **choi, int nlogs, double complex ***pauli);

/*
	Compute all the metrics for a given channel in the Choi matrix form.
	Inputs:
		double *metvals = double array of shape (nmetrics), which will contain the metric values after function execution.
		int nmetrics = number of metric values to be computed.
		char **metricsToCompute = array of strings containing the names of the metrics to be computed.
		double complex **choi = array of shape (4 x 4) containing the Choi matrix of the channel whose metrics need to be computed.
		double **ptm = Pauli transfer matrix.
		char *chname = name of the channel.
		struct constants_t *consts = Pointer to the constants structure, to access Pauli matrices.
*/
extern void ComputeMetrics(double *metvals, int nmetrics, char **metricsToCompute, double **ptm, char *chname, struct constants_t *consts);

#endif /* LOGMETRICS_H */