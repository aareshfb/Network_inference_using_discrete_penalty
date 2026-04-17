#!/bin/bash

# List of float values to send
t_VALUES=("10")
seed_VALUES=("1994" "1995" "1996" "1997" "1998")

# Loop through the nu_VALUES list and submit a job for each

for t_VAL in "${t_VALUES[@]}"
do
    for seed_VAL in "${seed_VALUES[@]}"
    do
        sbatch --job-name="jointGHS${seed_VAL}_${t_VAL}" exe_AUPRC.sh "$t_VAL" "$seed_VAL"
    done
done