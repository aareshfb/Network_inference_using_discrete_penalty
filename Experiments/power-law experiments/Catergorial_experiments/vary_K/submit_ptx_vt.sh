#!/bin/bash

# List of float values to send
t_VALUES=("3" "5" "10" "20" "50" "100")
seed_VALUES=("1337" "1338" "1339" "1340" "1341")
# nu_VALUES=("0.0001" "0.0005" "0.001" "0.005")

# Loop through the nu_VALUES list and submit a job for each

for t_VAL in "${t_VALUES[@]}"
do
    for seed_VAL in "${seed_VALUES[@]}"
    do
        sbatch --job-name="${seed_VAL}_${t_VAL//./}" exe_ptx_vt.sh "$t_VAL" "$seed_VAL"
    done
done