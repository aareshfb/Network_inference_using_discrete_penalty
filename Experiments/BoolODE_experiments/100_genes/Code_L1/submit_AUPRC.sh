#!/bin/bash

# List of float values to send
n_VALUES=("2500")


for n_VAL in "${n_VALUES[@]}"
do
    sbatch --job-name="L1${n_VAL//./}" exe_auprc.sh "$n_VAL"
    sleep 1
done
