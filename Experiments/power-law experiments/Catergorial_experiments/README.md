# Categorical Experiments (Power-Law)

This directory contains experiments for the **categorical (multi-population) setting** with **two categories**.

## Description

Each subdirectory corresponds to a different experimental setting:

- **`vary_K/`**: Experiments with varying number of populations \( K \)  
- **`vary_np/`**: Experiments with varying \( n/p \) ratio (number of samples relative to number of genes)  
- **`vary_perturb_strength/`**: Experiments with varying perturbation strength (`delta`) in the simulated data  

All experiments are conducted in the **two-category setting**.

## Running Experiments

To generate results for a specific setting, navigate to the desired subdirectory and run:

```bash
cd vary_perturb_strength    # or vary_K / vary_np
bash submit_ptx.sh 
```
## Notes

- These experiments evaluate the performance of the categorical algorithm under different data regimes.
- Scripts for running these experiments follow the same structure as those in `BoolODE_experiments/` and `scMTNI_experiments/`.
