#!/bin/bash

# List of float values to send
n_VALUES=("200" "1000" "2000")


for n_VAL in "${n_VALUES[@]}"
do
    sbatch --job-name="${n_VAL//./}" exe_auprc.sh "$n_VAL"
    sleep 1
done
