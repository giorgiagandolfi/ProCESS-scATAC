#!/bin/bash
#SBATCH --partition=smp
#SBATCH --job-name=simulate_chrAcc
#SBATCH --output=simulate_chrAcc.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1
#SBATCH --time=5:00:00
#SBATCH --mem-per-cpu=8GB


base="/data/rds/DMP/UCEC/GENEVOD/ggandolfi/Github/ProCESS-scATAC/scripts/simulate_fragments/peak_based/1_simulate_chromatin_profile/"

image="/data/rds/DMP/UCEC/GENEVOD/ggandolfi/Github/ProCESS-scATAC/scripts/simulate_fragments/peak_based/process_1.3_05.sif"


activity_list_file="/data/rds/DMP/UCEC/GENEVOD/ggandolfi/Github/ProCESS-scATAC/scripts/simulate_fragments/peak_based/input_data_P05/data/final_version/a_scores_pathway.rds"
peak_df_file="/data/rds/DMP/UCEC/GENEVOD/ggandolfi/Github/ProCESS-scATAC/scripts/simulate_fragments/peak_based/input_data_P05/data/final_version/peak_pathway_list.rds"


#### phylo forest path
sample_forest_file="/data/rds/DMP/UCEC/GENEVOD/ggandolfi/Github/ProCESS-scATAC/scripts/simulate_fragments/peak_based/0_process_simulations/sample_forest_atac_epigenome_1.3.5_pat05.sff"
phylo_forest_file="/data/rds/DMP/UCEC/GENEVOD/ggandolfi/Github/ProCESS-scATAC/scripts/simulate_fragments/peak_based/0_process_simulations/phylo_forest_atac_epigenome_1.3.5_pat05.sff"


#### outfiles
# out_cna_file="input_data_P05/data/df_peak_cna_final_gene.rds"
# out_peak_acc_file="input_data_P05/data/df_peak_final_gene.rds"
# 

srun singularity exec $image Rscript "$base/1_simulate_chromatin_accessibility_version_1.3.5.R" "$activity_list_file" "$peak_df_file" "$sample_forest_file" "$phylo_forest_file"
