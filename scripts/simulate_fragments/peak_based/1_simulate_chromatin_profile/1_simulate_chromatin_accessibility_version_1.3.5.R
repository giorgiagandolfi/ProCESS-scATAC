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
setwd("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-scATAC/scripts/simulate_fragments/peak_based")
source("utils.R")

sample_forest <- load_sample_forest("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-scATAC-old/scripts/simulate_fragments/peak_based/sample_forest_atac_epigenome_new_version_chloe_data.sff")
phylo_forest <- load_phylogenetic_forest("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-scATAC-old//scripts/simulate_fragments/peak_based/phylo_forest_atac_epigenetic_new_version_chloe_data.sff")


####### define the markov chain of open-closed chromatin
crc_peaks <- read_excel('/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-scATAC-old/scripts/public_data_processing/41586_2023_6682_MOESM4_ESM (1).xlsx',sheet = 1) %>%
  filter(Cancer=="CRC")
total_simulated_peaks <- nrow(crc_peaks)



activity <- list(
  'A' = list('E1'=c("P1"=1, "P2"=0.8, "P3"=0.9),
             'E2'=c("P1"=1, "P2"=0.9, "P3"=0.9),
             'E3'=c("P1"=1, "P2"=1, "P3"=1),
             'E4'=c("P1"=1, "P2"=0.9, "P3"=0.8)),
  'B' = list('E1'=c("P1"=1, "P2"=0.8, "P3"=0.9),
             'E2'=c("P1"=1, "P2"=0.9, "P3"=0.9),
             'E3'=c("P1"=1, "P2"=1, "P3"=1),
             'E4'=c("P1"=1, "P2"=0.9, "P3"=0.8)),
  'C' = list('E1'=c("P1"=1, "P2"=0.8, "P3"=0.9),
             'E2'=c("P1"=1, "P2"=0.9, "P3"=0.9),
             'E3'=c("P1"=1, "P2"=1, "P3"=1),
             'E4'=c("P1"=1, "P2"=0.9, "P3"=0.8)),
  'D' = list('E1'=c("P1"=1, "P2"=0.8, "P3"=0.9),
             'E2'=c("P1"=1, "P2"=0.9, "P3"=0.9),
             'E3'=c("P1"=1, "P2"=1, "P3"=1),
             'E4'=c("P1"=1, "P2"=0.9, "P3"=0.8))
)
get_epigenetic_activity<- function(activity,mutant,clone){
  programs <- activity[[mutant]][[clone]]
  return(programs)
}




# -----------------------------
# 3. Peaks + process assignment
# -----------------------------
n_peaks_P1 <- round(total_simulated_peaks/3,0)
n_peaks_P2 <- round(total_simulated_peaks/3,0)
n_peaks_P3 <- total_simulated_peaks-n_peaks_P1-n_peaks_P2

processes <- c(rep("P1",n_peaks_P1),rep("P2",n_peaks_P2),rep("P3",n_peaks_P3))

crc_peaks <- crc_peaks %>%
  mutate(process=processes) %>%
  separate(peak,into = c("chr","from","to"),sep = "-",remove = F,convert = T) %>%
  mutate(chr = str_remove(chr, "chr")) %>%
  mutate(peak_lenght=to-from)
sampled_lines <- sample(x = seq_along(1:nrow(crc_peaks)))
crc_peaks <- crc_peaks[sampled_lines,]

####### only for testing purposes
crc_peaks <- crc_peaks %>%
  add_row(
    Cancer = "CRC",
    peak = "chr7-120423815-120424315",
    chr = "7",
    from = 120423815,
    to = 120424315,
    pct.1 = 0.5,
    pct.2 = 0.2,
    p_val_adj = "1e-10",
    avg_log2FC = 1.5,
    process = "P1",
    peak_lenght = 500
  )
crc_peaks <- crc_peaks %>%
  select(peak,chr,from,to,peak_lenght,process)

selected_frags_dist_len = readRDS("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-scATAC-old/scripts/simulate_fragments/peak_based/selected_frags_dist_len.rds")
simulated_programs_name <- c("P1","P2","P3")
program_peak_list = lapply(simulated_programs_name,function(prog_id){
  crc_peaks %>% filter(process==prog_id)
})
names(program_peak_list)=simulated_programs_name

peak_access_labelling <- function(node) {
  cell_id_node = node$cell_id
  message(paste0("Cell: ",cell_id_node))

  cell_mutant <- node$mutant_name
  cell_phenotype <- node$epistate_name
  clone_programs <-get_epigenetic_activity(activity = activity,mutant = cell_mutant,clone=cell_phenotype)
  cell_activity_peaks <- list()
  for (p in names(clone_programs)){

    program_peaks <- program_peak_list[[p]]
    a_score = clone_programs[p]

    program_status = rbinom(nrow(program_peaks), size = 1, prob = a_score)

    cell_activity_peaks[[p]]<-program_peaks %>%
      mutate(
        status = program_status,
        cell_id = cell_id_node,
        mutant = cell_mutant,
        epistate = cell_phenotype
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

start=Sys.time()


sampled_cells = phylo_forest$get_nodes() %>% filter(!is.na(sample)) %>% pull(cell_id)
for (cell_id in sampled_cells[1]){
  node <- phylo_forest$get_node(cell_id)
  cell_peaks <- peak_access_labelling(node)
  genome <- node$get_genome()
  simulated_frags_list <- list()
  fasta_cell <- file(paste0("fasta_test_",cell_id,".txt"), open = "w")
  for (cp in 1:nrow(cell_peaks)){
    if (cell_peaks$status[cp]==1){
      cell_alleles = genome$get_alleles_covering_ref_region(cell_peaks$chr[cp],
                                                            cell_peaks$from[cp],
                                                            cell_peaks$peak_lenght[cp])
      simulated_frags_cell <- sample_fragments_for_peak_vec_allele(
        peak_id   = cell_peaks$peak[cp],
        peak_chr = cell_peaks$chr[cp],
        peak_from = as.numeric(cell_peaks$from[cp]),
        peak_to   = as.numeric(cell_peaks$to[cp]),
        fragment_len_dist = selected_frags_dist_len,
        available_alleles = cell_alleles,
        cell_id   = cell_id
      ) %>%  as.data.frame() %>%
        ungroup()
      simulated_frags_list[[cp]]<-simulated_frags_cell
      for (f in 1:nrow(simulated_frags_cell)){
        fragment_region <- genome$get_region_aligned_to_ref(simulated_frags_cell$fragment_chr[f],
                                                            simulated_frags_cell$fragment_allele[f]%>% as.numeric(),
                                                            simulated_frags_cell$fragment_start[f] %>% as.numeric(),
                                                            simulated_frags_cell$fragment_size[f] %>% as.numeric() )

        fragment = genome$get_fragment(fragment_region$chr,
                                       fragment_region$allele,
                                       fragment_region$from,
                                       fragment_region$length)
        fragment_seq = fragment$sequence
        fragment_id = paste0(">",cell_id,":",simulated_frags_cell$fragment_chr[f],":", simulated_frags_cell$fragment_start[f],":",simulated_frags_cell$fragment_allele[f])
        writeLines(fragment_id, fasta_cell)
        writeLines(fragment_seq, fasta_cell)

      }

    } else {
      message("Closed chromatin peak")
    }
  }
  simulated_frags_df = do.call("rbind",simulated_frags_list)
  close(fasta_cell)
}
end=Sys.time()
end-start
