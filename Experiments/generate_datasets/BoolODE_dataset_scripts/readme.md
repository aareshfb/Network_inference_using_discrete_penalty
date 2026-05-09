### Files

`simulate_GRNs.R` : main script where we set parameters for network evolution (detailed in Algorithm 2) and generate the ground-truth GRNS for all populations. Includes functionality for converting the adjacency matrix to boolean rules to input to BoolODE

`sparsify_counts.R`: Script to add (80%) sparsity to BoolODE simulated counts across all genes

`Initial-conditions.txt` : initial concentration of gene products for the BoolODE model. We always initialize the network with 1 for all regulators and 0 for target genes

`config.yaml` : BoolODE configuration file to set number of cells to simulate for each population, and number of time points to run simulations for
