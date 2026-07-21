# rm(list=ls())
library(ProCESS)
library(dplyr)
# library(ggplot2)
library(tidyverse)
library(readxl)
library(parallel)
# library(ComplexHeatmap)
# library(circlize)
# library(data.table)
#source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC/scripts/plot_muller.R")
#setwd("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-scATAC/scripts/simulate_fragments/peak_based/")
#source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-scATAC/scripts/simulate_fragments/peak_based/utils.R")
setwd("/data/rds/DMP/UCEC/GENEVOD/ggandolfi/Github/ProCESS-scATAC/scripts/simulate_fragments/peak_based/1_simulate_chromatin_profile")
source("../utils.R")

args <- commandArgs(trailingOnly = TRUE)
#### acticity and gene files
activity_list_file <- (args[1])
peak_df_file <- (args[2])
#### phylo forest path
sample_forest_file <- args[3]
phylo_forest_file <- args[4]


### read info
sample_forest <- load_sample_forest(sample_forest_file)
phylo_forest <- load_phylogenetic_forest(phylo_forest_file)
activity <- readRDS(activity_list_file)
peaks <-  readRDS(peak_df_file)



selected_frags_dist_len = readRDS("/data/rds/DMP/UCEC/GENEVOD/ggandolfi/Github/ProCESS-scATAC/scripts/simulate_fragments/peak_based/selected_frags_dist_len.rds")


peak_access_labelling <- function(node) {
  cell_id_node = node$cell_id
  message(paste0("Cell: ",cell_id_node))

  cell_mutant <- node$mutant_name
  cell_epistate <- node$epistate_name
  clone_programs <-get_epigenetic_activity(activity = activity,epistate=cell_epistate)
  cell_activity_peaks <- list()
  for (p in names(clone_programs)){

    program_peaks <- peaks[[p]]
    a_score = clone_programs[p]

    program_status = rbinom(nrow(program_peaks), size = 1, prob = a_score)

    cell_activity_peaks[[p]]<-program_peaks %>%
      mutate(
        status = program_status,
        cell_id = cell_id_node,
        mutant = cell_mutant,
        epistate = cell_epistate
      )

    cell_active_peaks <- do.call("rbind",cell_activity_peaks)
    return(cell_active_peaks)
  }
}



start=Sys.time()
node_tour <- get_node_tour(phylo_forest, only_leaves = T,with_genomes = T)
end=Sys.time()
end-start


###node =phylo_forest$get_node()
### node$

# start=Sys.time()
# parallel::detectCores()
# 
sampled_cells = phylo_forest$get_nodes() %>% filter(!is.na(sample)) %>% pull(cell_id)
# for (cell_id in sampled_cells[1:5]){
#   node <- phylo_forest$get_node(cell_id)
#   cell_peaks <- peak_access_labelling(node)
#   genome <- node$get_genome()
#   simulated_frags_list <- list()
#   fasta_cell <- file(paste0("fasta_test_",cell_id,".txt"), open = "w")
#   for (cp in 1:nrow(cell_peaks)){
#     if (cell_peaks$status[cp]==1){
#       cell_alleles = genome$get_alleles_covering_ref_region(cell_peaks$chr[cp],
#                                                             cell_peaks$from[cp],
#                                                             cell_peaks$peak_lenght[cp])
#       simulated_frags_cell <- sample_fragments_for_peak_vec_allele(
#         peak_id   = cell_peaks$peak[cp],
#         peak_chr = cell_peaks$chr[cp],
#         peak_from = as.numeric(cell_peaks$from[cp]),
#         peak_to   = as.numeric(cell_peaks$to[cp]),
#         fragment_len_dist = selected_frags_dist_len,
#         available_alleles = cell_alleles,
#         cell_id   = cell_id
#       ) %>%  as.data.frame() %>%
#         ungroup()
#       simulated_frags_list[[cp]]<-simulated_frags_cell
#       for (f in 1:nrow(simulated_frags_cell)){
#         fragment_region <- genome$get_region_aligned_to_ref(simulated_frags_cell$fragment_chr[f],
#                                                             simulated_frags_cell$fragment_allele[f]%>% as.numeric(),
#                                                             simulated_frags_cell$fragment_start[f] %>% as.numeric(),
#                                                             simulated_frags_cell$fragment_size[f] %>% as.numeric() )
# 
#         fragment = genome$get_fragment(fragment_region$chr,
#                                        fragment_region$allele,
#                                        fragment_region$from,
#                                        fragment_region$length)
#         fragment_seq = fragment$sequence
#         fragment_id = paste0(">",cell_id,":",fragment_region$chr,":", fragment_region$from,":",fragment_region$allele)
#         writeLines(fragment_id, fasta_cell)
#         writeLines(fragment_seq, fasta_cell)
# 
#       }
# 
#     } else {
#       message("Closed chromatin peak")
#     }
#   }
#   simulated_frags_df = do.call("rbind",simulated_frags_list)
#   saveRDS(object = simulated_frags_df,file = paste0('frags_',cell_id,".rds"))
#   close(fasta_cell)
# }

library(parallel)
library(dplyr)

# number of CPUs allocated by SLURM
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
if (is.na(ncores)) ncores <- 2  # fallback for RStudio/local runs

sampled_cells <- phylo_forest$get_nodes() %>%
  filter(!is.na(sample)) %>%
  pull(cell_id)

start <- Sys.time()

mclapply(sampled_cells, function(cell_id) {
  
  fasta_file <- paste0("cell_", cell_id, ".fasta")
  rds_file <- paste0("frags_cell_", cell_id, ".rds")
  
  # Skip already completed cells
  if (file.exists(fasta_file) && file.exists(rds_file)) {
    message("Skipping ", cell_id, " (files already exist)")
    return(cell_id)
  }
  
  node <- phylo_forest$get_node(cell_id)
  cell_peaks <- peak_access_labelling(node)
  genome <- node$get_genome()
  
  simulated_frags_list <- list()
  
  fasta_cell <- file(
    paste0("cell_", cell_id, ".fasta"),
    open = "w"
  )
  
  for (cp in 1:nrow(cell_peaks)) {
    
    if (cell_peaks$status[cp] == 1) {
      
      cell_alleles <- genome$get_alleles_covering_ref_region(
        cell_peaks$chr[cp],
        cell_peaks$from[cp],
        cell_peaks$peak_lenght[cp]
      )
      
      simulated_frags_cell <- sample_fragments_for_peak_vec_allele(
        peak_id = cell_peaks$peak[cp],
        peak_chr = cell_peaks$chr[cp],
        peak_from = as.numeric(cell_peaks$from[cp]),
        peak_to = as.numeric(cell_peaks$to[cp]),
        fragment_len_dist = selected_frags_dist_len,
        available_alleles = cell_alleles,
        cell_id = cell_id
      ) %>%
        as.data.frame() %>%
        ungroup()
      
      simulated_frags_list[[cp]] <- simulated_frags_cell
      
      for (f in 1:nrow(simulated_frags_cell)) {
        
        fragment_region <- genome$get_region_aligned_to_ref(
          simulated_frags_cell$fragment_chr[f],
          as.numeric(simulated_frags_cell$fragment_allele[f]),
          as.numeric(simulated_frags_cell$fragment_start[f]),
          as.numeric(simulated_frags_cell$fragment_size[f])
        )
        
        fragment <- genome$get_fragment(
          fragment_region$chr,
          fragment_region$allele,
          fragment_region$from,
          fragment_region$length
        )
        
        fragment_seq <- fragment$sequence
        
        fragment_id <- paste0(
          ">",
          cell_id,
          ":",
          fragment_region$chr,
          ":",
          fragment_region$from,
          ":",
          fragment_region$allele
        )
        
        writeLines(fragment_id, fasta_cell)
        writeLines(fragment_seq, fasta_cell)
      }
      
    } else {
      message("Closed chromatin peak")
    }
  }
  
  simulated_frags_df <- do.call(
    "rbind",
    simulated_frags_list
  )
  
  saveRDS(
    object = simulated_frags_df,
    file = paste0("frags_cell_", cell_id, ".rds")
  )
  
  close(fasta_cell)
  
  return(cell_id)
  
}, mc.cores = ncores)


end=Sys.time()
end-start
