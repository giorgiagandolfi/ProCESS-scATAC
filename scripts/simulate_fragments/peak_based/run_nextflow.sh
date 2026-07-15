#!/bin/bash
#SBATCH --job-name=nf_master
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --partition=master-worker
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4000
#SBATCH --time=72:00:00

job_id=$SLURM_JOB_ID

# ---------------------------------------------------------------------------
# Run the comprehensive sequenza pipeline
#   BAMs ??? seqz ??? first fit ??? refit at alternative solutions
#
# Usage:
#   sbatch run_pipeline.sh
#   # or locally:
#   bash run_pipeline.sh
#
# Edit MANIFEST, OUTDIR, REF_GENOME and GC_GENOME before submitting.
# ---------------------------------------------------------------------------

PIPELINE_DIR="/data/rds/DMP/UCEC/GENEVOD/ggandolfi/Github/ProCESS-scATAC/nf-atac_reads_simulator"
CONFIG="/data/rds/DMP/UCEC/GENEVOD/ggandolfi/Github/ProCESS-scATAC/nf-atac_reads_simulator/alma_giorgia.config"
OUTDIR="/data/rds/DMP/UCEC/GENEVOD/ggandolfi/Github/ProCESS-scATAC/scripts/simulate_fragments/peak_based/test_nextflow_pipeline"
INPUT="/data/rds/DMP/UCEC/GENEVOD/ggandolfi/Github/ProCESS-scATAC/scripts/simulate_fragments/peak_based/samplesheet.csv"
REF_GENOME="/data/rds/DMP/UCEC/GENEVOD/ggandolfi/Github/ProCESS-scATAC/scripts/simulate_fragments/peak_based/process_references_RACES/GRCh38/reference_tmp.fasta"
BOWTIE="/data/rds/DMP/UCEC/GENEVOD/ggandolfi/Github/ProCESS-scATAC/scripts/simulate_fragments/peak_based/bowtie2/reference_tmp"
# Load environment
source ~/.bashrc
conda activate my_nextflow

# Run the pipeline
nextflow run "${PIPELINE_DIR}/main.nf" \
    --input "${INPUT}" \
    -profile alma \
    -c "${CONFIG}" \
    --outdir "${OUTDIR}" \
    --fasta "${REF_GENOME}" \
    --bowtie2_index "${BOWTIE}"
