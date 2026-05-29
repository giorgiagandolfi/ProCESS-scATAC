#!/bin/bash
#SBATCH --partition=compute
#SBATCH --job-name=merge_bam
#SBATCH --output=merge_bam_%A_%a.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1
#SBATCH --time=3:00:00
#SBATCH --mem-per-cpu=4GB

# -----------------------
# Paths / container
# -----------------------
IMAGE1="/data/scratch/shared/SINGULARITY-DOWNLOAD/nextflow/.singularity/depot.galaxyproject.org-singularity-samtools-1.23--h96c455f_0.img"
# extract cell ID from filename
FILES=($(ls fragments_cells/cell_*.bam | sort))
OUT_BAM="fragments_cells/merged.bam"

srun singularity exec "$IMAGE1" \
    samtools merge "$OUT_BAM" "${FILES[@]}" 

srun singularity exec "$IMAGE1" \
    samtools index "$OUT_BAM"
