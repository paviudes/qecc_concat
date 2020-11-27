/*
	To compile this file, use the compiler options and link commands in file:///opt/intel/documentation_2019/en/mkl/common/mkl_link_line_advisor.htm .
	Do the following commands to compile this file.
		source /opt/intel/compilers_and_libraries_2019/mac/bin/compilervars.sh intel64
		gcc -m64 -I${MKLROOT}/include -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_rt -lpthread -lm -ldl mt19937/mt19937ar.c printfuns.c checks.c -o checks.o
*/
#include <stdlib.h>
#include <stdio.h>
#include <time.h> // Only for testing purposes
#include <string.h> // Only for testing purposes
#include <math.h>
#include <complex.h>
#include "constants.h"
#include "reps.h"
#include "mt19937/mt19937ar.h" // Only for testing purposes
#include "printfuns.h" // Only for testing purposes
#include "linalg.h"
#include "checks.h"

int IsDiagonal(long double **matrix, int size){
	// Check if the input matrix is diagonal or not.
	const long double atol = 10E-15;
	int i, j;
	for (i = 0; i < size; i ++)
		for (j = 0; j < size; j ++)
			if (i != j)
				if (fabsl(matrix[i][j]) > atol)
					return 0;
	return 1;
}

int IsPositive(double complex **choi, double atol){
	// Test if an input complex 4x4 matrix is completely positive.
	// A completely positive matrix has only non-negative eigenvalues.
	double complex *eigvals = malloc(sizeof(double complex) * 4);
	DiagonalizeD(choi, 4, eigvals, 0, NULL);
	int i, ispos = 1;
	for (i = 0; i < 4; i ++){
		if (cimag(eigvals[i]) > atol)
			ispos = 0;
		if (creal(eigvals[i]) < -1 * atol)
			ispos = 0;
	}
	if (ispos == 0){
		PrintComplexArray2D(choi, "Choi", 4, 4);
		PrintComplexArray1D(eigvals, "eigenvalues", 4);
	}
	// PrintComplexArray1D(eigvals, "eigenvalues", 4);
	// Free memory
	free(eigvals);
	return ispos;
}

int IsHermitian(double complex **choi, double atol){
	// Check is a complex 4x4 matrix is Hermitian.
	// For a Hermitian matrix A, we have: A[i][j] = (A[j][i])^*.
	int i, j;
	for (i = 0; i < 4; i ++){
		for (j = 0; j < 4; j ++){
			if (cabs(choi[i][j] - conj(choi[j][i])) > atol){
				return 0;
			}
		}
	}
	return 1;
}

int IsTraceOne(double complex **choi, double atol){
	// Check if the trace of a complex 4x4 matrix is 1.
	int i;
	double complex trace = 0 + 0 * I;
	for (i = 0; i < 4; i ++)
		trace = trace + choi[i][i];
	// printf("trace = %g + i %g.\n", creal(trace), cimag(trace));
	if (fabs(cimag(trace)) > atol){
		printf("trace = %.15f + i %.15f.\n", creal(trace), cimag(trace));
		return 0;
	}
	if (fabs(creal(trace)  -  1.0) > atol){
		printf("trace = %.15f + i %.15f.\n", creal(trace), cimag(trace));
		return 0;
	}
	return 1;
}

int IsState(double complex **choi, double atol){
	// Check is a 4 x 4 matrix is a valid density matrix.
	// It must be Hermitian, have trace 1 and completely positive.
	int isstate = 0, ishermitian = 0, ispositive = 0, istrace1 = 0;
	ishermitian = IsHermitian(choi, atol);
	if (ishermitian == 0)
		printf("Not Hermitian.\n");
	ispositive = IsPositive(choi, atol);
	if (ispositive == 0)
		printf("Not Positive.\n");
	// istrace1 = IsTraceOne(choi, atol);
	istrace1 = 1; // ONLY FOR DEBUGGING
	if (istrace1 == 0)
		printf("Not Unit trace.\n");
	isstate = ishermitian * istrace1 * ispositive;
	return isstate;
}

int _IsChannel(double **ptm, struct constants_t *consts, double atol){
	// Check is a 4 x 4 matrix is a valid density matrix.
	// We will convert it to a Choi matrix and test if the result is a density matrix.
	double complex **choi = NULL;
	int r, c, nlogs = 4;
	choi = malloc(sizeof(double complex *) * nlogs);
	for (r = 0; r < nlogs; r ++){
		choi[r] = malloc(sizeof(double complex) * nlogs);
		for (c = 0; c < nlogs; c ++){
			choi[r][c] = 0 + 0 * I;
		}
	}
	ProcessToChoi(ptm, choi, nlogs, consts->pauli);
	// PrintLongDoubleArray2D(ptm, "PTM", nlogs, nlogs);
	// PrintComplexArray2D(choi, "Choi", nlogs, nlogs);
	int ischan = IsState(choi, atol);
	if (ischan == 0)
		PrintDoubleArray2D(ptm, "PTM", nlogs, nlogs);
	for (r = 0; r < nlogs; r ++)
		free(choi[r]);
	free(choi);
	return ischan;
}

int IsChannel(long double **ptm, struct constants_t *consts, double atol){
	// Check is a 4 x 4 matrix is a valid density matrix.
	// We will convert it to a Choi matrix and test if the result is a density matrix.
	double **ptmd = malloc(sizeof(double *)*4);
	int i,j;
	for (i=0; i<4 ; i++){
		ptmd[i] = malloc(sizeof(double)*4);
		for (j=0; j<4 ; j++)
			ptmd[i][j] = (double) ptm[i][j];
	}
	int ischan = _IsChannel(ptmd, consts, atol);
	for (i=0; i<4 ; i++)
		free(ptmd[i]);
	free(ptmd);
	return ischan;
}


int IsPDF(long double *dist, int size){
	// Check if a given list of numbers is a normalized PDF.
	const long double atol = 10E-8;
	int i;
	long double norm = 0;
	for (i = 0; i < size; i ++){
		norm = norm + dist[i];
		if (dist[i] < 0){
			return 0;
		}
	}
	if (fabsl(norm) - 1 > atol){
		return 0;
	}
	return 1;
}

/*
int main(int argc, char const *argv[])
{
	init_genrand(time(NULL));
	int i, j, k;
	double complex **mat = malloc(sizeof(double complex) * 4);
	for (i = 0; i < 4; i ++){
		mat[i] = malloc(sizeof(double complex) * 4);
		for (j = 0; j < 4; j ++)
			mat[i][j] = genrand_real3() + genrand_real3() * I;
	}
	if (strncmp(argv[1], "IsPositive", 10) == 0){
		printf("Function: IsPositive.\n");
		PrintComplexArray2D(mat, "Matrix", 4, 4);
		int ispos = IsPositive(mat);
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
			ispos = IsPositive(matP);
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
		int isherm = IsHermitian(mat);
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
			isherm = IsHermitian(matH);
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
		int istrone = IsTraceOne(mat);
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
			istrone = IsTraceOne(matN);
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
		int isstate = IsState(mat);
		if (isstate == 0)
			printf("is not a state.\n");
		else
			printf("is a state.\n");
	}
	if (strncmp(argv[1], "IsPDF", 5) == 0){
		double sum = 0;
		double *dist = malloc(sizeof(double) * 100);
		for (i = 0; i < 100; i ++){
			dist[i] = genrand_real3();
			sum += dist[i];
		}
		PrintDoubleArray1D(dist, "Un-normalized distribution", 100);
		int ispdf = IsPDF(dist, 100);
		printf("has ispdf = %d.\n", ispdf);
		for (i = 0; i < 100; i ++)
			dist[i] = dist[i]/sum;
		PrintDoubleArray1D(dist, "And after normalization", 100);
		ispdf = IsPDF(dist, 100);
		printf("it has ispdf = %d.\n", ispdf);
		// Free memory
		free(dist);
	}
	// Free memory
	for (i = 0; i < 4; i ++)
		free(mat[i]);
	free(mat);
}
*/
