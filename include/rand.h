#ifndef RAND_H
#define RAND_H

// Compute an integer random number in a range.
extern int RandomRangeInt(int from, int to);

// Populate an integer array with random entries in a range.
extern void RandomArray1D(int from, int to, int size, int *arr);

// Shuffle the entries of the sequence.
extern void ShuffleInt(int *arr, int size, int nshuffles);

// Compute a random binary sequence of a given weight.
// First define a binary sequence of a given weight w.
// Shuffle the entries of the sequence.
extern void RandomBinarySequence(int nbits, int weight, int *arr);

#endif /* RAND_H */