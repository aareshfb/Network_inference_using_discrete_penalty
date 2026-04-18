# Categorical Experiment Results

This directory contains the results and post-processing scripts for categorical experiments under the power-law setting.


## Description

- **`1337/`, `1338/`, ..., `1341/`**:  
  Results from individual runs (different random seeds)

- **`Cat_results*.csv`**:  
  Results for the categorical method across runs  

- **`Direct_results*.csv`**:  
  Results for the direct (baseline) method  

- **`*_averaged.csv`**:  
  Averaged results across multiple runs  

- **`average_data.py`**:  
  Script to aggregate results across runs  

- **`generate_plots.py`**:  
  Script to generate plots from processed results  

- **`Plots/`**:  
  Contains generated figures used in the paper  

## Usage

To compute averaged results:

```bash
python3 average_data.py
```
