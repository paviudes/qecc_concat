#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <complex.h>
#include "printfuns.h"

void PrintComplexArray1D(complex double *array, char *name, int nrows){
	// Print a complex 2D array.
	int r;
	printf("-----\n");
	printf("Array name: %s\n", name);
	for (r = 0; r < nrows; r ++)
		printf("    %g + i %g", creal(array[r]), cimag(array[r]));
	printf("\n");
	printf("-----\n");
}

void PrintComplexArray2D(complex double **array, char *name, int nrows, int ncols){
	// Print a complex 2D array.
	int r, c;
	printf("-----\n");
	printf("Array name: %s\n", name);
	for (r = 0; r < nrows; r ++){
		for (c = 0; c < ncols; c ++)
			printf("    %g + i %g", creal(array[r][c]), cimag(array[r][c]));
		printf("\n");
	}
	printf("-----\n");
}


void PrintDoubleArrayDiag(double **array, char *name, int nrows){
	// Print the diagonal of a double 2D array.
	int r;
	printf("-----\n");
	printf("Diagonal array name: %s\n", name);
	for (r = 0; r < nrows; r ++){
		printf("    %.12f", (array[r][r]));
	}
	printf("\n");
	printf("-----\n");
}

void PrintDoubleArray2D(double **array, char *name, int nrows, int ncols){
	// Print a double 2D array.
	int r, c;
	printf("-----\n");
	printf("Array name: %s\n", name);
	for (r = 0; r < nrows; r ++){
		for (c = 0; c < ncols; c ++)
			printf("    %.16f", (array[r][c]));
		printf("\n");
	}
	printf("-----\n");
}

void CopyDoubleArray2D(double **from, double **to, int nrows, int ncols){
	// Copy the elements of one array into another.
	int r, c;
	// PrintDoubleArray2D(from, "From", nrows, ncols);
	for (r = 0; r < nrows; r ++){
		for (c = 0; c < ncols; c ++){
			to[r][c] = from[r][c];
		}
	}
	// PrintDoubleArray2D(to, "To", nrows, ncols);
}


void PrintDoubleArray1D(double *array, char *name, int nrows){
	// Print a double 1D array.
	int r;
	printf("-----\n");
	printf("Array name: %s\n", name);
	for (r = 0; r < nrows; r ++)
		printf("    %.16f", (array[r]));
	printf("\n");
	printf("-----\n");
}

void PrintIntArray2D(int **array, char *name, int nrows, int ncols){
	// Print a int 2D array.
	int r, c;
	printf("-----\n");
	printf("Array name: %s\n", name);
	for (r = 0; r < nrows; r ++){
		for (c = 0; c < ncols; c ++)
			printf("    %d", (array[r][c]));
		printf("\n");
	}
	printf("-----\n");
}

void PrintIntArray1D(int *array, char *name, int nrows){
	// Print a int 1D array.
	int r;
	printf("-----\n");
	printf("Array name: %s\n", name);
	for (r = 0; r < nrows; r ++)
		printf("    %d", (array[r]));
	printf("\n");
	printf("-----\n");
}

void PrintLongArray2D(long **array, char *name, int nrows, int ncols){
	// Print a long 2D array.
	int r, c;
	printf("-----\n");
	printf("Array name: %s\n", name);
	for (r = 0; r < nrows; r ++){
		for (c = 0; c < ncols; c ++)
			printf("    %ld", (array[r][c]));
		printf("\n");
	}
	printf("-----\n");
}

void PrintLongArray1D(long *array, char *name, int nrows){
	// Print a long 1D array.
	int r;
	printf("-----\n");
	printf("Array name: %s\n", name);
	for (r = 0; r < nrows; r ++)
		printf("    %ld", (array[r]));
	printf("\n");
	printf("-----\n");
}

void LoadDoubleArray1D(double *array, char *fname, int size){
	// Print a double 1D array.
	int i;
	if( access( fname, F_OK ) != -1 ){
		FILE *afp = fopen(fname, "r");
		for (i = 0; i < size; i ++)
			fscanf(afp, "%lf", &(array[i]));
		fclose(afp);
	}
	else
		printf("File %s does not exist.\n", fname);
}

void LoadIntArray1D(int *array, char *fname, int size){
	// Print a double 1D array.
	int i;
	if( access( fname, F_OK ) != -1 )
	{
		FILE *afp = fopen(fname, "r");
		for (i = 0; i < size; i ++)
			fscanf(afp, "%d", &(array[i]));
		fclose(afp);
	}
	else
		printf("File %s does not exist.\n", fname);
}