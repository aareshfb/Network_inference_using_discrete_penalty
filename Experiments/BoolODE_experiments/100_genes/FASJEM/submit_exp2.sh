#!/bin/bash

# List of float values to send
n_VALUES=("2500")

# Loop through the nu_VALUES list and submit a job for each

for n_VAL in "${n_VALUES[@]}"
do
    sbatch --job-name="FASJEM${n_VAL}" exe_auprc.sh "$n_VAL"
done