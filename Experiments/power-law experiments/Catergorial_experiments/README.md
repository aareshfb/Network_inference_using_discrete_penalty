# Categorical Experiments (Power-Law)

This directory contains experiments for the **categorical (multi-population) setting** with **two categories**.

## Description

Each subdirectory corresponds to a different experimental setting:

- **`vary_K/`**: Experiments with varying number of populations \( K \)  
- **`vary_np/`**: Experiments with varying \( n/p \) ratio (number of samples relative to number of genes)  
- **`vary_perturb_strength/`**: Experiments with varying perturbation strength (`delta`) in the simulated data  

All experiments are conducted in the **two-category setting**.

## Notes

- These experiments evaluate the performance of the categorical algorithm under different data regimes.
- Detailed descriptions of the implementation and files in this directory are provided in `vary_perturb_strength/README.md`. The same structure applies to the other settings (`vary_np/`,`vary_K/`).
- Scripts for running these experiments follow the same structure as those in `BoolODE_experiments/` and `scMTNI_experiments/`.
