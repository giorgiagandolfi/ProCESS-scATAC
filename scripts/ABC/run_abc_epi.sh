#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --job-name=process_abc_epi
#SBATCH --output=process_abc_epi.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --time=6:00:00
#SBATCH --mem-per-cpu=10GB


basedir='/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project'
base="${basedir}/ProCESS-scATAC/scripts/ABC/"
image="${basedir}/singularity_images/process_1_3_09_abc.sif"

module load singularity


singularity exec --bind /orfeo:/orfeo --no-home $image Rscript "$base/abc_single_clone_with_epi.R"
#srun singularity exec $image Rscript "$base/save_cna_cell.R"
