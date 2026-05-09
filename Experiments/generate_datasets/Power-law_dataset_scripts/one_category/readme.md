###Files

`simulation_scaleFree.R` : The script to generate random MST through k populations, and power-law distributed modular GRNs for each population. User can set number of genes (p), number of modules (M), number of populations (k), number of modules perturbed at each branch of population-graph (s) and the magnitude of perturbations (delta), number of cells to simulate for each population (n). The script generates the ground-truth GRNs and normally-distributed expression data with this precision-matrix for each population.

`simulation_helpers.R` : helper functions to run the main script.
