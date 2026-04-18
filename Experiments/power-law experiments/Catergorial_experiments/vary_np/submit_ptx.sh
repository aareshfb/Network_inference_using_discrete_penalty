#!/bin/bash

# List of float values to send
np_VALUES=("1" "5" "10" "15" "20" "25" "30")
seed_VALUES=("10195" "2338" "339" "40" "51341")
ptx_VALUES=("1.0")
######### nu_VALUES=("0.0001" "0.0005" "0.001" "0.005")

# Loop through the nu_VALUES list and submit a job for each

for np_VAL in "${np_VALUES[@]}"
do
    for seed_VAL in "${seed_VALUES[@]}"
    do
        for ptx_VAL in "${ptx_VALUES[@]}"
        do
            sbatch --job-name="${seed_VAL}_${np_VAL//./}" exe_ptx_vnp.sh "$np_VAL" "$seed_VAL" "$ptx_VAL"
        done
    done
done
