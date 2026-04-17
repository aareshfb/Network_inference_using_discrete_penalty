# BoolODE Experiments

This directory contains datasets and scripts for reproducing the BoolODE-based experiments corresponding to **Sim-1**, **Sim-2**, and **Sim-3** in the paper.

## Description

Each subdirectory corresponds to a simulation setting with a different number of genes:

- **Sim-1 (`65_genes/`)**: \( p = 65 \) genes  
- **Sim-2 (`100_genes/`)**: \( p = 100 \) genes  
- **Sim-3 (`200_genes/`)**: \( p = 200 \) genes  

Within each directory:
- Gene expression datasets generated using the **BoolODE** framework are provided.
- Algorithm-specific subdirectories contain scripts for:
  - Data loading  
  - Running the methods  
  - Computing evaluation metrics 
  - Bash scripts (`submit_AUPRC.sh`) specifying hyperparameters and launching experiments  

## Running Experiments

To reproduce results, navigate to a specific simulation setting and run the corresponding script in an algorithm subdirectory:

```bash
cd 65_genes/Code_L0
bash submit_AUPRC.sh
```
## Note

- Additional implementation details are provided in the `65_genes/README.md` file. The same setup applies to the other simulation settings (`100_genes/` and `200_genes/`).
