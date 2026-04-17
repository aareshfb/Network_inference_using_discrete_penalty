#!/bin/sh
#SBATCH --job-name=F"$1"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:15:00
#SBATCH --partition=standard
#SBATCH --licenses=gurobi@slurmdb:1
#SBATCH --account=fattahi0

module load Rtidyverse/4.4.0
module load RStudio/2024.04.1
 

Rscript run_SILGGM_AUPRC.R "$1" > code_print.txt