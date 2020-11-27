#ifndef CHECKS_H
#define CHECKS_H

#include <complex.h>

// Check is a matrix is a diagonal matrix or not.
extern int IsDiagonal(long double **matrix, int size);

// Check if a given list of numbers is a normalized PDF.
extern int IsPDF(long double *dist, int size);

// Check if an input complex 4x4 matrix is completely positive.
extern int IsPositive(double complex **choi, double atol);

// Check if the trace of a complex 4x4 matrix is 1.
extern int IsTraceOne(double complex **choi, double atol);

// Check is a complex 4x4 matrix is Hermitian.
extern int IsHermitian(double complex **choi, double atol);

// Check is a 4 x 4 Choi matrix is a valid density matrix.
extern int IsState(double complex **choi, double atol);

// Check is a 4 x 4 Pauli transfer matrix is a valid channel.
extern int _IsChannel(double **ptm, struct constants_t *consts, double atol);
extern int IsChannel(long double **ptm, struct constants_t *consts, double atol);

#endif /* CHECKS_H */
