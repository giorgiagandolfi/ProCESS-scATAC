#!/bin/bash
#SBATCH --partition=smp
#SBATCH --job-name=process_epigenetic
#SBATCH --output=process_epigenetic.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=10GB

#image='/data/rds/DMP/UCEC/GENEVOD/ggandolfi/process_on_the_fly_v2.sif'
#base='/data/rds/DMP/UCEC/GENEVOD/ggandolfi'

base="/data/rds/DMP/UCEC/GENEVOD/ggandolfi/Github/ProCESS-scATAC/scripts/simulate_fragments/peak_based/0_process_simulations"
image="/data/rds/DMP/UCEC/GENEVOD/ggandolfi/Github/ProCESS-scATAC/scripts/simulate_fragments/peak_based/process_label_tour_v0.sif"
image="/data/rds/DMP/UCEC/GENEVOD/ggandolfi/Github/ProCESS-scATAC/scripts/simulate_fragments/peak_based/process_1.3_v1.sif"


srun singularity exec $image Rscript "$base/simulate_process_epigenetic_pat05.R"

