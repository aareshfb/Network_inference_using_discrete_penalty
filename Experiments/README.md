# Description of Folders

#### `BoolODE_experiments/` 
Contains scripts, datasets and experiment results based on **BoolODE-generated datasets**.  
These experiments evaluate performance on biologically motivated synthetic data with known ground-truth networks.

#### `Datasets/`
Contains the scripts used to generate the datasets used in the experiments.  

This includes:
- Power-law datasets   
- BoolODE datasets  

Each dataset is formatted as gene expression matrices used as input to the methods.

#### `power-law experiments/`
Contains scripts for experiments on **synthetic datasets generated using power-law network structures**.  

#### `scMTNI_experiments/`
Contains scripts for experiments on datasets from the **scMTNI benchmark**.  
These provide a more realistic evaluation setting for gene regulatory network inference.

## Notes

- Each experiment folder includes scripts for:
  - Data loading  
  - Running methods  
  - Evaluating performance (e.g., AUPRC, F1 score)  

- For synthetic datasets, results are typically averaged over multiple runs to account for randomness.

- All experiments can be reproduced by running the corresponding scripts in each directory.

## References

If you use the scMTNI datasets or benchmarking framework, please cite:

- Zhang, S. et al. *Inference of cell type-specific gene regulatory networks on cell lineages using multi-task learning*. Nature Communications, 2023.

- scMTNI GitHub repository: https://github.com/Roy-lab/scMTNI
