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
python3 code_category_vperturb.py "$nu_VAL" "True" > code_print.txt
```

- **`nu_VAL`**: parameter controlling the model
- **`True`**: flag for experiment settings (see script for details)
- **Note:** The hyperparameters `lambda` (sparsity) and `gamma` (similarity) are specified within `run_boolode_AUPRC.py`.


## Running Experiments

To generate results, run:

```bash
bash submit_ptx.sh
```
This script launches experiments with predefined parameter settings.

##Notes
- This setup is consistent with other experiment directories (e.g., BoolODE and scMTNI).

