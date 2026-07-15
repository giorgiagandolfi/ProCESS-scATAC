#!/bin/bash
#SBATCH --partition=compute
#SBATCH --job-name=call_peaks
#SBATCH --output=call_peaks_%A_%a.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1
#SBATCH --time=3:00:00
#SBATCH --mem-per-cpu=4GB

# -----------------------
# Paths / container
# -----------------------
IMAGE1="/data/rds/DMP/UCEC/GENEVOD/ggandolfi/Github/ProCESS-scATAC/scripts/simulate_fragments/peak_based/my_singularity_image/sinto%3A0.10.1--pyhdfd78af_0"

IN_BAM="fragments_cells/merged.bam"
OUT_BED="fragments_cells/called_peaks.bed"


srun singularity exec "$IMAGE1" \
    Genrich -j -t "$IN_BAM" \
    -o "$OUT_BED" -v