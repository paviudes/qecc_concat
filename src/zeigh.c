void ZEigH(double complex **mat, int dim, double complex *eigvals, int iseigvecs, double complex **eigvecs){
	/*
		Compute the eigenvalues and the right-eigenvectors of a complex Hermitian matrix.
		We will use the LAPACK routine zheev to compute the eigenvalues. The LAPACK function is defined as follows.
		extern void zgeev(char* jobvl, // Should left eigenvectors be computed? 'V' for yes and 'N' for no.
						  char* jobvr, // Should right eigenvectors be computed? 'V' for yes and 'N' for no.
						  int* n, // The order of the matrix.
						  dcomplex* a, // complex array containing the n x n complex matrix in vectorized form.
						  int* lda, // The leading dimension of the matrix.
						  dcomplex* w, // Array of the computed eigenvalues. Complex array of size n.
						  dcomplex* vl, // The left eigenvectors of A, stored as columns of this n x n matrix. The order of the eigenvectors is the same as the order of the eigenvalues.
						  int* lvdl, // The leading dimension of the array vl.
						  dcomplex* vr, // The right eigenvectors of A, stored as columns of this n x n matrix. The order of the eigenvectors is the same as the order of the eigenvalues.
						  int* lvdr, // The leading dimension of the array vr.
						  dcomplex* work, // A scratch workspace for the algorithm.
						  int* lwork, // Size of the array: work. If it is -1, then running the algorithm doesn't compute the eigenvalues but assigns the optimal size of the work array to lwork.
						  double* rwork, // double precision array of size 2 * n
						  int* info); // result. 0 if successful, -i if i-th argument has illegal value and +i if 1 to i eigen values of w have not converged.

		There is a C wrapper to the LAPACK routine:
		https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/mkl_lapack_examples/lapacke_zgeev_row.c.htm
	*/

	MKL_INT n = dim, lda = dim, info;
	MKL_Complex16 w[dim];
	MKL_Complex16 a[dim * dim];
	int i, j;
	for (i = 0; i < dim; i ++){
		for (j = 0; j < dim; j ++){
			a[i * dim + j].real = (creal(mat[i][j]));
			a[i * dim + j].imag = (cimag(mat[i][j]));
		}
	}
	for (i = 0; i < dim; i ++){
		for (j = 0; j < i; j ++){
			a[i * dim + j].real = 0;
			a[i * dim + j].imag = 0;
		}
	}

	info = LAPACKE_zheev(LAPACK_ROW_MAJOR, 'V', 'L', n, a, lda, w);
	if (info > 0)
		printf("Eigenvalues %d to %d did not converge properly.\n", info, dim);
	else if (info < 0)
		printf("Error in the the %d-th input parameter.\n", -1 * info);
	else{
		for (i = 0; i < dim; i ++){
			eigvals[i] = w[i].real + w[i].imag * I;
			if (iseigvecs == 1){
				for (j = 0; j < dim; j ++)
					eigvecs[i][j] = a[i * dim + j].real + a[i * dim + j].imag * I;
			}
		}
	}
}