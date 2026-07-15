#!/bin/bash
#SBATCH --partition=compute
#SBATCH --job-name=bowtie2_align
#SBATCH --output=bowtie2_align_%A_%a.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1
#SBATCH --time=3:00:00
#SBATCH --mem-per-cpu=4GB
#SBATCH --array=501-1143%200   # replace N=1143
# -----------------------
# Paths / container
# -----------------------
IMAGE1="/data/scratch/shared/SINGULARITY-DOWNLOAD/nextflow/.singularity/depot.galaxyproject.org-singularity-bowtie2-2.4.4--py39hbb4e92a_0.img"
IMAGE2="/data/scratch/shared/SINGULARITY-DOWNLOAD/nextflow/.singularity/depot.galaxyproject.org-singularity-samtools-1.23--h96c455f_0.img"
BOWTIE_REF="bowtie2/reference_tmp"
INDIR="fragments_cells_big"
FILES=($(ls ${INDIR}/cell_*.bed | sort))
echo $BED
BED=${FILES[$SLURM_ARRAY_TASK_ID]}

# extract cell ID from filename
CELL=$(basename "$BED" .bed)

CELL_FA="${INDIR}/fastqs/${CELL}.fa"
OUT_FASTQ="${INDIR}/fastqs//${CELL}.pe.fastq"

R1_FASTQ="${INDIR}/fastqs/${CELL}.R1.fastq"
R2_FASTQ="${INDIR}/fastqs/${CELL}.R2.fastq"

OUTDIR="${INDIR}/aligned_bam"
mkdir -p ${OUTDIR}
OUT_SAM="${OUTDIR}/${CELL}.sam"
OUT_CB_SAM="${OUTDIR}/${CELL}.cb.sam"
OUT_BAM="${OUTDIR}/${CELL}.bam"


echo "Processing $CELL"
echo "BED: $BED"
echo "OUT: $OUT"mac

srun singularity exec "$IMAGE1" \
   bowtie2   --very-sensitive   -X 2000 \
   -x "$BOWTIE_REF" \
   -1 "$R1_FASTQ" \
   -2 "$R2_FASTQ" \
   --rg-id "$CELL" \
   --rg "SM:$CELL" \
   -S "$OUT_SAM"

awk -v CB="$CELL" '
BEGIN { OFS="\t" }
/^@/ { print; next }
{ print $0, "CB:Z:" CB }
' "$OUT_SAM" > "$OUT_CB_SAM"

srun singularity exec "$IMAGE2" \
  samtools sort -o "$OUT_BAM" "$OUT_CB_SAM"

rm "$OUT_SAM"
rm "$OUT_CB_SAM"
srun singularity exec "$IMAGE2" \
  samtools index "$OUT_BAM"

