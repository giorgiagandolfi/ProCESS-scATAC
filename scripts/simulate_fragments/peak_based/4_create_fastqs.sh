#!/bin/bash
#SBATCH --partition=compute
#SBATCH --job-name=create_fastq
#SBATCH --output=create_fastq_%A_%a.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1
#SBATCH --time=3:00:00
#SBATCH --mem-per-cpu=4GB
#SBATCH --array=0-5 # replace N
# -----------------------
# Paths / container
# -----------------------
IMAGE1="/data/rds/DMP/UCEC/GENEVOD/ggandolfi/Github/ProCESS-scATAC/scripts/simulate_fragments/peak_based/my_singularity_image/art_modern-1.3.2.sif"
IMAGE2="/data/rds/DMP/UCEC/GENEVOD/ggandolfi/Github/ProCESS-scATAC/scripts/simulate_fragments/peak_based/my_singularity_image/bbmap%3A39.81--h9b5c0a0_1"

INDIR="fragments_cells_big"
FILES=($(ls ${INDIR}/cell_*.bed | sort))
echo $BED
BED=${FILES[$SLURM_ARRAY_TASK_ID]}

# extract cell ID from filename
CELL=$(basename "$BED" .bed)

CELL_FA="${INDIR}/reference_fasta/${CELL}.fa"

OUTDIR="${INDIR}/fastqs"
mkdir -p ${OUTDIR}
OUT_FASTQ="${OUTDIR}/${CELL}.pe.fastq"

R1_FASTQ="${OUTDIR}/${CELL}.R1.fastq"
R2_FASTQ="${OUTDIR}/${CELL}.R2.fastq"
echo "Processing $CELL"
echo "BED: $BED"
echo "OUT: $OUT"

srun singularity exec "$IMAGE1" \
    art_modern   --mode template \
    --lc pe \
    --i-file "$CELL_FA" \
    --o-fastq "$OUT_FASTQ" --i-fcov 2   --read_len 50 --min_qual 5
    
srun singularity exec "$IMAGE2" \
    reformat.sh in="$OUT_FASTQ" \
    out1="$R1_FASTQ" \
    out2="$R2_FASTQ"

rm "$OUT_FASTQ" 
