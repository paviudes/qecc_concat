void SetFullProcessMatrix(struct qecc_t *qecc, struct simul_t *sim, double *process, int isPauli){
	// In cases where the full process matrix needs to be explicity specified,
	// load the structure attribute with the specified matrix. In the case of
	// Pauli channels, the process matrix input must be a 1 x 4^n, each column
	// specifying a unique Pauli error probability, where the errors are ordered
	// as SLT.
	
}