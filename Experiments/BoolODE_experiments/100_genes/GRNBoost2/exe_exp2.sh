#!/bin/sh
#SBATCH --job-name=s"$3"_n"$2"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:15:00
#SBATCH --partition=standard
#SBATCH --licenses=gurobi@slurmdb:1
#SBATCH --account=fattahi0

module load python3.9-anaconda/2021.11 
conda init bash
source ~/.bashrc

conda activate conda1-env 


python3 Run_GRNBoost2.py > code_print.txt