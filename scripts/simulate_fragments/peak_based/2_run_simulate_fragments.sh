#!/bin/bash
#SBATCH --partition=compute
#SBATCH --job-name=simulate_fragments
#SBATCH --output=simulate_fragments_%A_%a.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --time=3:00:00
#SBATCH --mem-per-cpu=4GB
#SBATCH --array=0-5 # replace N
# -----------------------
# Paths / container
# -----------------------
IMAGE="/data/rds/DMP/UCEC/GENEVOD/ggandolfi/Github/ProCESS-scATAC/scripts/simulate_fragments/peak_based/process_label_tour_v0.sif"
base="/data/rds/DMP/UCEC/GENEVOD/ggandolfi/Github/ProCESS-scATAC/scripts/simulate_fragments/peak_based"
CELL_IDX=$((SLURM_ARRAY_TASK_ID + 1))
srun singularity exec $IMAGE Rscript "$base/2_simulate_fragments.R" "$CELL_IDX"