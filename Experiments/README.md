## Description of Folders

### `BoolODE_experiments/`
Contains scripts and experiments based on **BoolODE-generated datasets**.  
These experiments evaluate performance on biologically motivated synthetic data with known ground-truth networks.

### `Datasets/`
Contains all datasets used in the experiments.  

This includes:
- Synthetic datasets (e.g., power-law simulations)  
- BoolODE datasets  
- scMTNI benchmark datasets  

Each dataset is formatted as gene expression matrices used as input to the methods.

### `power-law experiments/`
Contains scripts for experiments on **synthetic datasets generated using power-law network structures**.  
These experiments are used to evaluate performance under controlled simulation settings.

### `scMTNI_experiments/`
Contains scripts for experiments on datasets from the **scMTNI benchmark**.  
These provide a more realistic evaluation setting for gene regulatory network inference.

## Notes

- Each experiment folder includes scripts for:
  - Data loading  
  - Running methods  
  - Evaluating performance (e.g., AUPRC, F1 score)  

- For synthetic datasets, results are typically averaged over multiple runs to account for randomness.

- All experiments can be reproduced by running the corresponding scripts in each directory.
