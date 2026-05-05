#!/bin/bash
#SBATCH --partition=GENOA
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem 80gb
#SBATCH --time=20:00:00
#SBATCH --output=ProCESS.out
#SBATCH --error=ProCESS.err

#module load singularity
image="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/process_1.1.0.sif"
# change with your path to the simulate_tissue.R and simulate_mutation.R scripts
base="/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC/scripts"
#
#singularity exec --bind /orfeo:/orfeo --no-home $image Rscript $base/simulate_tissue.R
#singularity exec --bind /orfeo:/orfeo --no-home $image Rscript $base/simulate_mutation.R

### test new ProCESS version - on-the-fly_mutations
module load R/4.4.1
Rscript $base/simulate_seq.R
