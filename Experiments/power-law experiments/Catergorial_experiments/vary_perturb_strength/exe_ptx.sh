#!/bin/sh
#SBATCH --job-name=s"$3"_n"$2"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:30:00
#SBATCH --partition=standard
#SBATCH --licenses=gurobi@slurmdb:1
#SBATCH --account=fattahi0

module load python3.9-anaconda/2021.11 gurobi/10.0.2
#conda activate conda-env-MIQP1  


ptx_VAL=$1
seed_VAL=$2

python3 code_category_vperturb.py "$ptx_VAL" 5 250 20 "$seed_VAL" > code_print.txt