# Simulation Setup (BoolODE - Sim-1)

This directory contains the datasets and implementation of all methods used in **Sim-1** (BoolODE-based experiments with \( p = 65 \) genes).

## Description

- **`Code_L0/`**: Implementation of **ELEM-0** (proposed method)  
- **`Code_L1/`**: Implementation of **ELEM-1**  
- **`FASJEM/`, `FGL/`, `GGL/`, `GRNBoost2/`, `SILGGM/`, `jointGHS/`**: Implementations of competing methods  
- **`Expression_data/`**: Gene expression datasets generated using BoolODE  
- **`Ground_truth/`**: Ground-truth networks used for evaluation  

Each method directory contains:
- Scripts for running the method  
- Evaluation scripts (e.g., AUPRC computation)  
- A `submit_AUPRC.sh` script specifying hyperparameters and launching experiments  

## Running Experiments

To run a specific method, navigate to its directory and execute:

```bash
cd Code_L0
bash submit_AUPRC.sh #Simulation Setup (BoolODE - Sim-1)
```

## Notes

- Detailed implementation and usage instructions are provided in `Code_L0/README.md` (ELEM-0).
- The same workflow and directory structure apply to the other methods (`Code_L1/`, `FASJEM/`, `FGL/`, `GGL/`, `GRNBoost2/`, `SILGGM/`, `jointGHS/`) unless otherwise noted.
