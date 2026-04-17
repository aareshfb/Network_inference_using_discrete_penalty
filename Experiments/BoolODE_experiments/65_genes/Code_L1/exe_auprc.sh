#!/bin/sh
#SBATCH --job-name=s"$3"_n"$2"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=1:30:00
#SBATCH --partition=standard
#SBATCH --licenses=gurobi@slurmdb:1
#SBATCH --account=fattahi0

# Positional arguments
n=$1

module load matlab
matlab -batch "n=$n; run('Run_L1_auprc_mask.m')"