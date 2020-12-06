#ifndef UTILS_H
#define UTILS_H

// #include <complex.h>

/*
	Compute the order of magnitude of a number to a base.
	Given x, compute y such as x = O(b^y).
	This should be called for x < 1.
*/
extern int OrderOfMagnitude(long double x, int b);

// Check if a number x is above the cut-off 10^c.
extern int IsAboveCutoff(long double x, int c);

// Compute the sign of a number, i.e., k such that (-1)^k = x.
extern int Sign(double x);

/*
	Compute the combinatorial factor: n choose k.
	C(n, k) = [n * (n-1) * .... * (n-k+1)] / [k * (k-1) * .... * 1]
*/
extern long Comb(int n, int k);

// Compute the maximum of two numbers.
extern double Max(double a, double b);

// Compute the minimum of two numbers.
extern double Min(double a, double b);

// Compute the sum of numbers in an array.
extern int SumInt(int *arr, int size);

// Divide small numbers.
extern long double Divide(long double num, long double den);

// Compute the sum of numbers in an array.
extern double SumDouble(double *arr, int size);
extern long double SumLongDouble(long double *arr, int size);

// Normalize the array.
extern void Normalize(long double *arr, int size);

// Compute the parity of the integer's binary representation.
extern int BitParity(int num);

/*
	Compute dot product modulo 2, between two integer's binary representations.
	https://stackoverflow.com/questions/43300462/most-efficient-way-to-evaluate-a-binary-scalar-product-mod-2
*/
extern int BinaryDot(int a, int b);

#endif /* UTILS_H */
