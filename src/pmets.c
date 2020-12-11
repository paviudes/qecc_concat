void ComputeLevelZeroMetrics(struct simul_t *sim, int nqubits, int nlogs, struct constants_t *consts){
	// Compute the level-0 (physical) metrics.
	long double **ptm;
	double *chanmets;
	int i, j, q;

	// printf("nlogs = %d, nqubits = %d\n", nlogs, nqubits);

	if (sim->iscorr == 0){
		// ProcessToChoi((sim->logical)[0], nlogs, choi, consts->pauli);
		ComputeMetrics((sim->metricValues)[0], sim->nmetrics, sim->metricsToCompute, (sim->logical)[0], sim->chname, consts);
	}
	else if (sim->iscorr == 2){
		// Compute the level 0 metrics for each qubit's physical channel and average them.
		// We are assuming that the error metrics under with tensor products.
		// Extract the single qubit maps from the physical channels.
		// printf("Computing Level 0 metrics\n");
		// Fidelity is multiplicative whereas the other error-metrics are additive.
		// Infidelity is 1 - the product of fidelities
		for (i = 0; i < sim->nmetrics; i ++)
			if (strncmp((sim->metricsToCompute)[i], "infid", 5) == 0)
				(sim->metricValues)[0][i] = 1;

		// PrintDoubleArray1D((sim->metricValues)[0], "Level-0 metrics", sim->nmetrics);
		ptm = malloc(sizeof(long double *) * nlogs);
		for (i = 0; i < nlogs; i++)
			ptm[i] = malloc(sizeof(long double) * nlogs);

		chanmets = malloc(sizeof(double) * nqubits);
		for (q = 0; q < nqubits; q ++){
			// Extract the physical channel for the qubit q.
			for (i = 0; i < nlogs; i ++){
				for (j = 0; j < nlogs; j ++){
					ptm[i][j] = (long double) ((sim->physical)[q * nlogs * nlogs + i * nlogs + j]);
				}
			}
			// PrintDoubleArray2D(ptm, "Pauli transfer matrix", nlogs, nlogs);

			// Compute metrics for the q-th qubit's channel.
			for (i = 0; i < sim->nmetrics; i ++)
				chanmets[i] = 0;
			ComputeMetrics(chanmets, sim->nmetrics, sim->metricsToCompute, ptm, sim->chname, consts);
			// printf("qubit: %d\n", q);
			// PrintDoubleArray1D(chanmets, "Physical metrics", sim->nmetrics);
			// Add the metrics across qubits for all except fidelity.
			// For fidelity, we need to multiply the values.
			for (i = 0; i < sim->nmetrics; i ++){
				// printf("Metric: %s\n", sim->metricsToCompute[i]);
				if (strncmp(sim->metricsToCompute[i], "infid", 5) == 0)
					(sim->metricValues)[0][i] *= (1 - chanmets[i]);
				else
					(sim->metricValues)[0][i] += chanmets[i];
			}
		}
		// Infidelity is 1 - the product of fidelities
		for (i = 0; i < sim->nmetrics; i ++)
			if (strncmp(sim->metricsToCompute[i], "infid", 5) == 0)
				(sim->metricValues)[0][i] = 1 - (sim->metricValues)[0][i];
		// Free memory
		for (i = 0; i < nlogs; i++)
			free(ptm[i]);
		free(ptm);
		free(chanmets);
	}
	else{};
	// PrintDoubleArray1D((sim->metricValues)[0], "Level-0 metrics", sim->nmetrics);
}