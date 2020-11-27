#include <stdlib.h>
#include <stdio.h>
// #include <time.h> // Only for testing purposes
// #include <string.h> // Only for testing purposes
#include <math.h>
// #include <complex.h>
// #include "printfuns.h" // Only for testing purposes
#include "utils.h"

/*
	To compile this file, use the compiler options and link commands in file:///opt/intel/documentation_2019/en/mkl/common/mkl_link_line_advisor.htm .
	Do the following commands to compile this file.
		source /opt/intel/compilers_and_libraries_2019/mac/bin/compilervars.sh intel64
		gcc mt19937/mt19937ar.c printfuns.c -m64 -I${MKLROOT}/include -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_rt -lpthread -lm -ldl linalg.c -o linalg.o
*/

long Comb(int n, int k){
	// Compute the combinatorial factor: n choose k.
	// C(n, k) = [n * (n-1) * .... * (n-k+1)] / [k * (k-1) * .... * 1]
	if (k > n - k){
		k = n - k;
	}
	double fact = 1;
	int i;
	for (i = 0; i < k; ++i){
		fact *= ((double) (n - i))/((double) (k - i));
	}
	return (long) fact;
}

double Max(double a, double b){
	// Compute the maximum of two numbers.
	// https://www.geeksforgeeks.org/conditional-or-ternary-operator-in-c-c/
	return (a > b) ? a:b;
}

double Min(double a, double b){
	// Compute the minimum of two numbers.
	// https://www.geeksforgeeks.org/conditional-or-ternary-operator-in-c-c/
	return (a < b) ? a:b;
}

int Sign(double x){
	// Compute the sign of a number, i.e., k such that (-1)^k = x.
	const double atol = 1E-15;
	if (fabs(x) > atol){
		// printf("Sign of %.15f is %d.\n", x, (x > atol) ? 1:-1);
		if (x > atol)
			return 1;
		else
			return -1;
	}
	// printf("Sign of %.15f = 0.\n", x);
	return 0;
}

int SumInt(int *arr, int size){
	// Compute the sum of numbers in an array.
	int i, sum = 0;
	for (i = 0; i < size; i ++)
		sum += arr[i];
	return sum;
}

long double SumLongDouble(long double *arr, int size){
	// Compute the sum of numbers in an array.
	int i;
	long double sum = 0;
	for (i = 0; i < size; i ++)
		sum += arr[i];
	return sum;
}

double SumDouble(double *arr, int size){
	// Compute the sum of numbers in an array.
	int i;
	double sum = 0;
	for (i = 0; i < size; i ++)
		sum += arr[i];
	return sum;
}

void Normalize(long double *arr, int size){
	// Normalize the array.
	// https://stackoverflow.com/questions/18069269/normalizing-a-list-of-very-small-double-numbers-likelihoods
	int i;
	long double *logarr = malloc(sizeof(long double) * size);
	long double maxlog = log10l(arr[0]);
	for (i = 0; i < size; i ++){
		logarr[i] = log10l(arr[i]);
		if (maxlog <= logarr[i])
			maxlog = logarr[i];
	}
	for (i = 0; i < size; i ++)
		logarr[i] = logarr[i] - maxlog;
	for (i = 0; i < size; i ++)
		arr[i] = pow(10, logarr[i]);

	long double sum = SumLongDouble(arr, size);
	for (i = 0; i < size; i ++)
		arr[i] = arr[i]/sum;

	free(logarr);
}

int BitParity(int num){
	// Compute the parity of the integer's binary representation.
	int parity = 0;
	if (num > 0){
		int i, nbits = (int)(1 + (log(num)/log(2)));
		// printf("log(%d) = %g, nbits = %d.\n", num, log(num), nbits);
		int *seq = malloc(sizeof(int) * nbits);
		for (i = nbits - 1; i >= 0; i --){
			seq[i] = (int)(num/pow(2, i));
			num = num % (int)(pow(2, i));
		}
		parity = SumInt(seq, nbits) % 2;
		free(seq);
	}
	return parity;
}

long double Divide(long double num, long double den){
	/*
		Divide small numbers.
		To avoid precision issues, we will perform this division in three steps.
		1. If the numerator's first N digits are "0", return 0.
		2. If the numerator and the denominator are alike up to first N digits, return 1.
		3. Do the log division otherwise.
	*/
	// const int digits = 15;
	long double result = 0;
	long double *mantissas = malloc(sizeof(long double) * 2);
	int *exps = malloc(sizeof(int) * 2);
	long double sign = (long double) (Sign(num) * Sign(den));
	/*
	if (fabsl(num) <= powl(10, -1 * (digits + 1)))
		result = 0;
	else if (fabsl(fabsl(num) - fabsl(den)) <= powl(10, -1 * digits))
		result = sign;
	else{
		mantissas[0] = frexpl(fabsl(num), &(exps[0]));
		mantissas[1] = frexpl(fabsl(den), &(exps[1]));
		// printf("A = %Lf x 2^%d and B = %Lf x 2^%d\n", mantissas[0], exps[0], mantissas[1], exps[1]);
		result = sign * ldexpl(mantissas[0]/mantissas[1], exps[0] - exps[1]);
	}
	*/
	mantissas[0] = frexpl(fabsl(num), &(exps[0]));
	mantissas[1] = frexpl(fabsl(den), &(exps[1]));
	// printf("A = %Lf x 2^%d and B = %Lf x 2^%d\n", mantissas[0], exps[0], mantissas[1], exps[1]);
	result = sign * ldexpl(mantissas[0]/mantissas[1], exps[0] - exps[1]);
	free(exps);
	free(mantissas);
	return result;
}

int BinaryDot(int a, int b){
	// Compute dot product modulo 2, between two integer's binary representations.
	// https://stackoverflow.com/questions/43300462/most-efficient-way-to-evaluate-a-binary-scalar-product-mod-2
	return BitParity(a & b);
}
