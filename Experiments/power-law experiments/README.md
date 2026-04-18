# Power-Law Experiments

This directory contains experiments and results for datasets generated using **power-law network structures**.

## Description

- **`Categorical_experiments/`**:  
  Contains scripts and outputs for experiments in the **multi-category (categorical) setting**.

  This includes:
  - Scripts for running the methods  
  - Evaluation scripts (e.g., F1 score computation)  
  - Output files used in the paper  

- **`Single_category/`**:  
  Contains output files for experiments in the **single-category setting**. Scripts follow the same structure as in `scMTNI_experiments/` and `BoolODE_experiments/` and are not duplicated.
  
## Notes

- The power-law datasets are used to evaluate performance under controlled synthetic settings.
- For experiments involving randomness, results are averaged over 5 runs.
