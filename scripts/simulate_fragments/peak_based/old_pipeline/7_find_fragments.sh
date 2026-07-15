#!/bin/bash
#SBATCH --partition=compute
#SBATCH --job-name=find_fragments
#SBATCH --output=find_fragments_%A_%a.log
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
IMAGE2="/data/scratch/shared/SINGULARITY-DOWNLOAD/nextflow/.singularity/depot.galaxyproject.org-singularity-tabix-1.11--hdfd78af_0.img"

IN_BAM="fragments_cells_big/merged.bam"
OUT_BED="fragments_cells_big/called_fragments.bed"
OUT_SORTED_BED="fragments_cells_big/called_fragments.sorted.tsv"
OUT_SORTED_BED_GZ="fragments_cells_big/called_fragments.sorted.tsv.gz"
srun singularity exec "$IMAGE1" \
    sinto fragments --bam "$IN_BAM" \
    -f "$OUT_BED" -t CB --use_chrom ".*"

sort -k1,1 -k2,2n "$OUT_BED" > "$OUT_SORTED_BED"

srun singularity exec "$IMAGE2" \
  bgzip -@ 4 "$OUT_SORTED_BED"

srun singularity exec "$IMAGE2" \ 
  tabix -p bed "$OUT_SORTED_BED_GZ"

# clean up
rm "$OUT_BED"

