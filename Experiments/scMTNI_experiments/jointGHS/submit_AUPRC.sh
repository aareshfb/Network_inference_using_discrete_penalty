#!/bin/bash

# List of float values to send
n_VALUES=("200" "1000" "2000")

# Loop through the nu_VALUES list and submit a job for each

for n_VAL in "${n_VALUES[@]}"
do
    sbatch --job-name="scMTNI-jointGHS${n_VAL}" exe_auprc.sh "$n_VAL"
done