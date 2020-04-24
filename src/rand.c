#include <stdlib.h>
#include <stdio.h>
#include "mt19937/mt19937ar.h"

int RandomRangeInt(int from, int to)
{
	// Compute an integer random number in a range.
	unsigned long randint = genrand_int32();
	return (int)((unsigned long)from + randint % (unsigned long)(to - from));
}

void RandomArray1D(int from, int to, int size, int *arr)
{
	// Populate an integer array with random entries in a range.
	int i;
	for (i = 0; i < size; i++)
		arr[i] = RandomRangeInt(from, to);
}

void ShuffleInt(int *arr, int size, int nshuffles)
{
	// Shuffle the entries of the sequence.
	int i, from, to, temp;
	for (i = 1; i < nshuffles; i++)
	{
		from = RandomRangeInt(1, size);
		to = RandomRangeInt(1, size);
		printf("%d). From %d to %d.\n", i, from, to);
		temp = arr[from];
		arr[from] = arr[to];
		arr[to] = temp;
	}
}

void RandomBinarySequence(int nbits, int weight, int *arr)
{
	// Compute a random binary sequence of a given weight.
	// First define a binary sequence of a given weight w.
	// Shuffle the entries of the sequence.
	int i;
	for (i = 0; i < weight; i++)
		arr[i] = 1;
	for (i = weight; i < nbits; i++)
		arr[i] = 0;
	ShuffleInt(arr, nbits, nbits);
}