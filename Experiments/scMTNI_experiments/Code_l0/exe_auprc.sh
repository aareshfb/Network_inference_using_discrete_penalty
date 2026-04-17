#!/bin/sh
#SBATCH --job-name=l0_n"$1"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=3:30:00
#SBATCH --partition=standard
#SBATCH --licenses=gurobi@slurmdb:1
#SBATCH --account=fattahi0

module load python3.9-anaconda/2021.11 gurobi/10.0.2
#conda activate conda-env-MIQP1  


nu_VAL=$1
n=$2

python3 run_scMTNI_AUPRC.py "$n" "$nu_VAL" "True" > code_print.txt