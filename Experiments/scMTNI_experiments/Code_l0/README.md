# ELEM-0 (scMTNI datasets 1-3)

This directory contains the implementation and scripts for running the **ELEM-0** algorithm on scMTNI data (datasets 1-3).

## Description

- **`ELEM_0.py`**: Script for the ELEM-0 algorithm.  
- **`Para_algo2.py`, `Break_points_fun2.py`**: Scripts related to the underlying optimization algorithm.  
- **`funcs.py`**: Supporting utility functions used by ELEM-0.  
- **`run_scMTNI_AUPRC.py`**: Runs a single instance of the algorithm for given parameters.  
- **`get_AUPRC.py`**: Processes the output and computes the AUPRC metric.  
- **`submit_AUPRC.sh`**: Submits jobs and sets parameter values (e.g., backward map parameters).  
- **`exe_auprc.sh`**: Helper script that calls `run_scMTNI_AUPRC.py`.  
- **`Output/AUPRC/`**: Contains output files used for evaluation and reporting results.

## Running a Single Experiment

To run a single instance of the algorithm:

```bash
python3 run_scMTNI_AUPRC.py "$nu_VAL" "True" > code_print.txt
```

- **`nu_VAL`**: parameter controlling the model
- **`True`**: flag for experiment settings (see script for details)
- **Note:** The hyperparameters `lambda` (sparsity) and `gamma` (similarity) are specified within `run_boolode_AUPRC.py`.

## Running Full Experiments

To reproduce all experiments, run:

```bash
bash submit_AUPRC.sh
```

This script:
- Sets backward mapping paramter `nu`
- Launches multiple runs via exe_auprc.sh

## Notes
- After running experiments, use get_AUPRC.py to compute evaluation metrics.
- Output files are stored in Output/AUPRC/.
- This setup is consistent with other methods (`FASJEM`,`JGL`, etc) in the repository.



