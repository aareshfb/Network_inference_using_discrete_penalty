#!/bin/sh
#SBATCH --job-name=JGL_%a
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=01:15:00
#SBATCH --partition=standard
#SBATCH --licenses=gurobi@slurmdb:1
#SBATCH --account=fattahi0
#SBATCH --output=logs/JGL_%A_%a.out
#SBATCH --error=logs/JGL_%A_%a.err

mkdir -p logs

t_VAL=10
seed_VALUES=(1994 1995 1996 1997 1998)

# number of lambda points
N_LAMBDA=40

lambda_VALUES=($(awk -v n=$N_LAMBDA 'BEGIN{
    for(i=0;i<n;i++){
        x = -6*i/(n-1);
        printf("%.12g ", 10^x);
    }
}'))

task_id=${SLURM_ARRAY_TASK_ID}

seed_idx=$((task_id / N_LAMBDA))
lam_idx=$((task_id % N_LAMBDA))

seed_VAL=${seed_VALUES[$seed_idx]}
lam_VAL=${lambda_VALUES[$lam_idx]}

echo "task_id = $task_id"
echo "t       = $t_VAL"
echo "seed    = $seed_VAL"
echo "lambda  = $lam_VAL"

module load Rtidyverse/4.4.0
module load RStudio/2024.04.1

Rscript run_JGL_AUPRC.R "$t_VAL" "$seed_VAL" "$lam_VAL"