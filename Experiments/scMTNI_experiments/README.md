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
