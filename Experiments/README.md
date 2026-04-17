# Description of Folders

#### `BoolODE_experiments/` 
Contains scripts, datasets and experiment results based on **BoolODE-generated datasets**.  
These experiments evaluate performance on biologically motivated synthetic data with known ground-truth networks.

#### `Datasets/`
Contains the scripts used to generate the datasets used in the experiments.  
This includes:
- Power-law datasets   
- BoolODE datasets  

#### `power-law experiments/`
Contains results on **synthetic datasets generated using power-law network structures**.  

#### `scMTNI_experiments/`
Contains the data, results and scripts for experiments on the **scMTNI benchmark**.  


## Notes

- The experiment folder contains algorithm specific subdirectories with scripts for:
  - Bash scripts that specify hyperparameters and launch the experiments
  - Data loading  
  - Running methods  
  - Calculating the error metrics (eg., AUPRC and F1 score) 

- For synthetic datasets generated using power-law, results are typically averaged over 5 runs to account for randomness.

- BoolODE and scMTNI experiments can be reproduced by running the `submit_AUPRC.sh` scripts located within each algorithm-specific subdirectory (e.g., `Code_L0/`, `FGL/`, `GRNBoost2/`).

## References

If you use the scMTNI datasets or benchmarking framework, please cite:

- Zhang, S., Pyne, S., Pietrzak, S. et al. Inference of cell type-specific gene regulatory networks on cell lineages from single cell omic datasets. Nat Commun 14, 3064 (2023). https://doi.org/10.1038/s41467-023-38637-9

- scMTNI GitHub repository: https://github.com/Roy-lab/scMTNI
