#!/bin/bash

# List of float values to send
ptx_VALUES=("0.1" "0.2" "0.05" "0.5" "1.0" "2.0")
seed_VALUES=("1337" "1338" "1339" "1340" "1341")
# nu_VALUES=("0.0001" "0.0005")

# Loop through the nu_VALUES list and submit a job for each

for ptx_VAL in "${ptx_VALUES[@]}"
do
    for seed_VAL in "${seed_VALUES[@]}"
    do
        sbatch --job-name="${seed_VAL}_${ptx_VAL//./}" exe_ptx.sh "$ptx_VAL" "$seed_VAL"
    done
done