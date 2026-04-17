# scMTNI Experiments

This directory contains datasets and scripts for reproducing experiments based on the **scMTNI benchmark**.


## Description

- **`Code_L0/`**: Implementation of **ELEM-0** (proposed method)  
- **`Code_L1/`**: Implementation of **ELEM-1**  
- **Other method folders** (`FASJEM/`, `FGL/`, `GGL/`, `GRNBoost2/`, `SILGGM/`, `jointGHS/`): Implementations of competing methods  
- **`Simulated_GRNs/`**: Contains simulated network structures used in the benchmark  

Each method directory contains:
- Scripts for running the method  
- Evaluation scripts (e.g., AUPRC computation)  
- A `submit_AUPRC.sh` script specifying hyperparameters and launching experiments  
- Output files used in the paper
- For additional details on the files within each directory, refer to `Code_L0/README.md`. The organization and workflow are identical.

## Running Experiments

To reproduce results, navigate to a method directory and run:

```bash
cd Code_L0
bash submit_AUPRC.sh
```
Replace `Code_L0` with the desired method.

## Notes

- Additional data files present in this directory correspond to the scMTNI datasets (including ground-truth networks) and are used directly by the scripts.
- The directory structure is consistent across methods.

## References

If you use the scMTNI datasets or benchmarking framework, please cite:

- Zhang, S., Pyne, S., Pietrzak, S. et al. Inference of cell type-specific gene regulatory networks on cell lineages from single cell omic datasets. Nat Commun 14, 3064 (2023). https://doi.org/10.1038/s41467-023-38637-9

- scMTNI GitHub repository: https://github.com/Roy-lab/scMTNI
