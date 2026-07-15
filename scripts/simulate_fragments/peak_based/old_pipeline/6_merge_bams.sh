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
FILES=($(ls fragments_cells_big/aligned_bam/cell_*.bam | sort))
OUT_BAM="fragments_cells_big/merged900_cells.bam"
BAM_LIST="fragments_cells_big_bam_list.txt"

srun singularity exec "$IMAGE1" \
    samtools merge "$OUT_BAM" "${FILES[@]:0:900}" 
#srun singularity exec "$IMAGE1" \
#    samtools merge "$OUT_BAM" -b "$BAM_LIST"
srun singularity exec "$IMAGE1" \
    samtools index "$OUT_BAM"
