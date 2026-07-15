#!/bin/bash
#SBATCH --partition=smp
#SBATCH --job-name=simulate_chrAcc
#SBATCH --output=simulate_chrAcc.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1
#SBATCH --time=3:00:00
#SBATCH --mem-per-cpu=8GB


base="/data/rds/DMP/UCEC/GENEVOD/ggandolfi/Github/ProCESS-scATAC/scripts/simulate_fragments/peak_based/1_simulate_chromatin_profile/"

image="/data/rds/DMP/UCEC/GENEVOD/ggandolfi/Github/ProCESS-scATAC/scripts/simulate_fragments/peak_based/process_label_tour_v0.sif"


activity_list_file="input_data_P05/data/activity_list_gene_level_log1pscore.rds"
peak_df_file="input_data_P05/data/all_peaks_gene_df.rds"


#### phylo forest path
phylo_forest_file="phylo_forest_atac_epigenetic.sff"


#### outfiles
out_cna_file="input_data_P05/data/df_peak_cna_final_gene.rds"
out_peak_acc_file="input_data_P05/data/df_peak_final_gene.rds"


srun singularity exec $image Rscript "$base/1_simulate_chromatin_acc_new.R" "$activity_list_file" "$peak_df_file" "$phylo_forest_file" "$out_cna_file" "$out_peak_acc_file"

