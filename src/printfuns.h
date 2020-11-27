#ifndef PRINTFUNS_H
#define PRINTFUNS_H

#include <complex.h>

// Print a complex 1D array.
extern void PrintComplexArray1D(double complex *array, char *name, int nrows);
// Print a long complex 1D array.
extern void PrintLongComplexArray1D(long double complex *array, char *name, int nrows);

// Print a complex 2D array.
extern void PrintComplexArray2D(double complex **array, char *name, int nrows, int ncols);
// Print a long complex 2D array.
extern void PrintLongComplexArray2D(long double complex **array, char *name, int nrows, int ncols);

// Print a double 2D array.
extern void PrintDoubleArray2D(double **array, char *name, int nrows, int ncols);
// Print a long double 2D array.
extern void PrintLongDoubleArray2D(long double **array, char *name, int nrows, int ncols);

// Print a column of a double 2D array.
extern void PrintDoubleColumn(double **array, char *name, int nrows, int col);

// Print a row of a double 2D array.
extern void PrintDoubleRow(double **array, char *name, int ncols, int row);

// Print a double 1D array.
extern void PrintDoubleArray1D(double *array, char *name, int nrows);
// Print a long double 1D array.
extern void PrintLongDoubleArray1D(long double *array, char *name, int nrows);

// Print a int 2D array.
extern void PrintIntArray2D(int **array, char *name, int nrows, int ncols);

// Print a int 1D array.
extern void PrintIntArray1D(int *array, char *name, int nrows);

// Print a long 2D array.
extern void PrintLongArray2D(long **array, char *name, int nrows, int ncols);

// Print a long 1D array.
extern void PrintLongArray1D(long *array, char *name, int nrows);

// Load a double 1D array from a file.
extern void LoadDoubleArray1D(double *array, char *fname, int size);

// Load a int 1D array from a file.
extern void LoadIntArray1D(int *array, char *fname, int size);

// Print the diagonal of a double 2D array.
extern void PrintDoubleArrayDiag(double **array, char *name, int nrows);
// Print the diagonal of a long double 2D array.
extern void PrintLongDoubleArrayDiag(long double **array, char *name, int nrows);

// Copy the elements of one array into another.
extern void CopyDoubleArray2D(double **from, double **to, int nrows, int ncols);

#endif /* PRINTFUNS_H */
