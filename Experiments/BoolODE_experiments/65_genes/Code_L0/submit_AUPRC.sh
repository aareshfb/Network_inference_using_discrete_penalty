#!/bin/bash

# List of float values to send

#nu_VALUES=("0.00001" "0.00005" "0.0001" "0.0005" "0.001" "0.005" "0.01" "0.05")
# nu_VALUES=("0.005")
nu_VALUES=("0.0005")
# Loop through the nu_VALUES list and submit a job for each


for nu_VAL in "${nu_VALUES[@]}"
do
  sbatch --job-name="${nu_VAL//./}_n${n}" exe_auprc.sh "$nu_VAL" "$n"
  sleep 1
done
