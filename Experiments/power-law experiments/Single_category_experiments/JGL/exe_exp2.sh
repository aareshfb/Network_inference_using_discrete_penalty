#!/bin/sh
#SBATCH --job-name=JGL"$3"_n"$2"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:15:00
#SBATCH --partition=standard
#SBATCH --licenses=gurobi@slurmdb:1
#SBATCH --account=fattahi0

t_VAL=$1
seed_VAL=$2

module load Rtidyverse/4.4.0
module load RStudio/2024.04.1
#conda activate conda-env-MIQP1  

Rscript run_JGL.R "$t_VAL" 250 20 "$seed_VAL" > code_print.txt