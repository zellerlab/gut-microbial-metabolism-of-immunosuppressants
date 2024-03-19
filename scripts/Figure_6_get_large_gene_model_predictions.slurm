#!/bin/bash -e
#SBATCH -A zeller
#SBATCH --time=0-02:00:00
#SBATCH --qos=low
## Looks like we need more than 100 GB of ram for the analysis...
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
### The following number might need to be updated!
#SBATCH --array=1-47975
#SBATCH --error=/scratch/karcher/tmp2/%a.err
#SBATCH --output=/scratch/karcher/tmp2/%a.out
#SBATCH --job-name=large_gene_models_full_models_all_drugs

#cat make_tmps/tmp8.18.sh_MODULE | sed "s/.*; //" |  head -n $SLURM_ARRAY_TASK_ID | tail -1 | bash
cd /g/scb/zeller/karcher/gut-microbial-metabolism-of-critical-dose-immunosuppressants/scripts
cat large_gene_model_lines.iv |  head -n $SLURM_ARRAY_TASK_ID | tail -1 | bash
