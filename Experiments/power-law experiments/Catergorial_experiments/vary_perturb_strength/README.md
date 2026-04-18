# Categorical Experiment (Power-Law Setting)

This directory contains scripts and outputs for running the categorical (two-category) experiments under a specific simulation setting.

## Description

- **`code_category_vperturb.py`**: Main script implementing the categorical algorithm  
- **`Para_algo2.py`, `Break_points_fun2.py`**: Scripts for the underlying optimization procedure  
- **`funcs.py`**: Supporting utility functions  
- **`load_results2.py`**: Processes outputs for evaluation  
- **`submit_ptx.sh`**: Submits jobs and sets experiment parameters  
- **`exe_ptx.sh`**: Helper script that executes individual runs  
- **`Output/`**: Contains output files used for evaluation and reporting results  

## Running a Single Experiment

To run a single instance of the algorithm:

```bash
python3 code_category_vperturb.py "$ptx_VAL" "$K_VAL" "$p_VAL" "$np_VAL" "$seed_VAL" > code_print.txt
```

- **`ptx_VAL`**: parameter controlling the perturb strenght `delta
- **`K_VAL`**: number of poulations
- **`p_VAL`**: number of genes
- **`np_VAL`**: number of samples to genes
- **Note:** The hyperparameters `lambda` (sparsity) and `gamma` (similarity) are specified within `run_boolode_AUPRC.py`.


## Running Experiments

To generate results, run:

```bash
bash submit_ptx.sh
```
This script launches experiments with predefined parameter settings.

- Random seeds used to generate the data and the perturb strength paramter are set in `submit_ptx.sh`
- **Limitation:** This functionality is not implemented in this dicrectory. To use it, the corresponding data must be added to the appropriate folder, and the relevant lines in `funcs.py` and `code_category_vperturb.py` must be modified accordingly.
  - Update the file path to data in `funcs.py` (line ~20) before running.
  - Update the file path to data in `code_category_vperturb.py` (line ~34,35,66) before running.
  - Data can be generated using the scrpits stored at `Experiments/Datsets/Power-law based datasets`
  
## Notes
- This setup is consistent with other experiment directories (e.g., BoolODE and scMTNI).

