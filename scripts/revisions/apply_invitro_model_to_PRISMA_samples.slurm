#!/bin/bash -e
#SBATCH -A zeller
#SBATCH --time=0-01:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
### The following number might need to be updated!
#SBATCH --array=1-57
#SBATCH --error=/scratch/karcher/tmp3/%a.err
#SBATCH --output=/scratch/karcher/tmp3/%a.out
#SBATCH --job-name=apply_large_model
cat iv.sh | head -n $SLURM_ARRAY_TASK_ID | tail -1 | bash
