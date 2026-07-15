#!/bin/bash
#SBATCH --partition=compute
#SBATCH --job-name=get_cell_refs
#SBATCH --output=get_cell_refs_%A_%a.log
#SBATCH --error=get_cell_refs_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --time=3:00:00
#SBATCH --mem-per-cpu=4GB
#SBATCH --array=401-1149%100 # replace N
# -----------------------
# Paths / container
# -----------------------
IMAGE="/data/rds/DMP/UCEC/GENEVOD/ggandolfi/Github/ProCESS-scATAC/scripts/simulate_fragments/peak_based/my_singularity_image/bedtools%3A2.31.1--hf5e1c6e_2"
FASTA="/data/rds/DMP/UCEC/GENEVOD/ggandolfi/Github/ProCESS-scATAC/scripts/simulate_fragments/peak_based/process_references_RACES/GRCh38/reference_tmp.fasta"

# get sorted list of all cell BED files
INDIR="fragments_cells_big"
FILES=($(ls ${INDIR}/cell_*.bed | sort))
BED=${FILES[$SLURM_ARRAY_TASK_ID]}

# extract cell ID from filename
CELL=$(basename "$BED" .bed)

OUTDIR="${INDIR}/reference_fasta"
mkdir -p ${OUTDIR}
OUT="${OUTDIR}/${CELL}.fa"

echo "Processing $CELL"
echo "BED: $BED"
echo "OUT: $OUT"

srun singularity exec "$IMAGE" \
    bedtools getfasta -name \
    -fi "$FASTA" \
    -bed "$BED" \
    -fo "$OUT"
