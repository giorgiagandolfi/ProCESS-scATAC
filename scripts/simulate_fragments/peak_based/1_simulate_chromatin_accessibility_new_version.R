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
source("utils.R")

sample_forest <- load_sample_forest("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-scATAC-old/scripts/simulate_fragments/peak_based/sample_forest_atac_epigenome_new_version.sff")
phylo_forest <- load_phylogenetic_forest("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-scATAC-old/scripts/simulate_fragments/peak_based/phylo_forest_atac_epigenetic_new_version.sff")


####### define the markov chain of open-closed chromatin
crc_peaks <- read_excel('/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-scATAC-old/scripts/public_data_processing/41586_2023_6682_MOESM4_ESM (1).xlsx',sheet = 1) %>%
  filter(Cancer=="CRC") %>%
  filter(grepl("chr7", peak))

total_simulated_peaks <- nrow(crc_peaks)

# activity <- list(
#   'A' = list('E1'=c("P1"=0.9, "P2"=0.8, "P3"=0.2),
#              'E2'=c("P1"=0.1, "P2"=0.3, "P3"=0.3),
#              'E3'=c("P1"=0.4, "P2"=0.6, "P3"=0.5),
#              'E4'=c("P1"=0.1, "P2"=0.2, "P3"=0.8))
# )


activity <- list(
  'A' = list('E1'=c("P1"=1, "P2"=0.8, "P3"=0.9),
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
sampled_lines <- sample(x = seq_along(1:nrow(crc_peaks)),size = 10)
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

cna_data_sim <- phylo_forest$get_cell_allelic_fragmentation()
cna_data_sim = cna_data_sim %>% 
  mutate(chr=paste0('chr',chr))



selected_frags_dist_len = readRDS("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-scATAC-old/scripts/simulate_fragments/peak_based/selected_frags_dist_len.rds")


random_labelling <- function(label,node) {
  cell_id_node = node$cell_id
  message(paste0("Cell: ",cell_id_node))
  cell_mutant = sample_forest$get_nodes() %>% 
    filter(cell_id==cell_id_node) %>% 
    pull(mutant)
  cell_phenotype = sample_forest$get_nodes() %>% 
    filter(cell_id==cell_id_node) %>% 
    pull(epistate)
  cell_cna = phylo_forest$get_sampled_cell_CNAs() %>% 
    filter(cell_id==cell_id_node) %>% 
    dplyr::rename(from=begin,
                  to=end) 
  
  cell_segments = phylo_forest$get_cell_allelic_fragmentation() %>% 
    filter(cell_id==cell_id_node) %>% 
    dplyr::rename(from=begin,
                  to=end) %>% 
    left_join(cell_cna)
  
  if (nrow(cell_segments)==0){
    message("Not sampled cell")
  } else {
    # cell_clone <- paste0(cell_mutant,cell_phenotype)
    clone_programs <-get_epigenetic_activity(activity = activity,mutant = cell_mutant,clone=cell_phenotype)
    cell_activity_peaks <- list()
    for (p in names(clone_programs)){
      program_peaks <- crc_peaks %>% filter(process==p)
      a_score = clone_programs[p]
      
      program_status = rbinom(nrow(program_peaks), size = 1, prob = a_score)
      
      peaks_dt <- as.data.table(program_peaks)
      # dplyr::rename(peak_chr=chr,
      #               peak_from=from,
      #               peak_to=to)
      cna_dt   <- as.data.table(cell_segments) %>% 
        dplyr::rename(seg_chr=chr,
                      seg_from=from,
                      seg_to=to)
      setkey(peaks_dt, chr, from, to)
      setkey(cna_dt, seg_chr, seg_from, seg_to)
      res =foverlaps(
        peaks_dt,
        cna_dt,
        by.x = c("chr", "from", "to"),
        by.y = c("seg_chr", "seg_from", "seg_to"),
        type = "any",   # overlap allowed
        nomatch = NA
      )
      df_peak_with_cna <- as.data.frame(res)
      
      cell_activity_peaks[[p]]<-program_peaks %>% 
        left_join(df_peak_with_cna) %>% 
        mutate(
          status = program_status,
          cell_id = cell_id_node,
          mutant = cell_mutant,
          epistate = cell_phenotype
        ) %>% 
        rowwise() %>% 
        mutate(
          available_alleles = list(
            case_when(
              is.na(type) ~ list(c(0, 1)),
              type == "A" ~ list(c(0, 1, allele)),
              type == "D" ~ list(src.allele),
              TRUE ~ list(numeric(0))
            )[[1]]
          )
        ) %>% 
        ungroup()
      
      
    }
    cell_active_peaks <- do.call("rbind",cell_activity_peaks) %>% 
      filter(status==1) %>% 
      mutate(tot_cn=minor+major)
    chromosomes <- cell_active_peaks %>% 
      pull(chr) %>% 
      unique()
    
    chromosomes_frags <- mclapply(chromosomes, mc.cores = detectCores() - 1,
                                  function(chrom){
                                    message(paste0("chromosome: ",chrom))  
                                    
                                    
                                    fragm_df_cell <- sample_fragments_for_peak_vec(
                                      peak_id   = cell_active_peaks$peak,
                                      peak_chr=str_replace(chrom,pattern = "chr",replacement = ""),
                                      peak_from = as.numeric(cell_active_peaks$from),
                                      peak_to   = as.numeric(cell_active_peaks$to),
                                      fragment_len_dist = selected_frags_dist_len,
                                      tot_cn    = as.numeric(cell_active_peaks$tot_cn),
                                      cell_id   = cell_id_node
                                    ) %>%  as.data.frame()
                                  }) %>% bind_rows()
    df_final = left_join(chromosomes_frags,cell_active_peaks)
    return(df_final)    
  }

  
}
# genome <- get_genome_tour(phylo_forest, only_leaves = TRUE)
label_tour <- get_label_tour(phylo_forest, random_labelling,
                             only_leaves = T,with_genomes = T)

list_peaks_final <- list()
i=1


while (!label_tour$done) {
  cell_id = label_tour$value$cell_id
  cell_frags = label_tour$value$label %>% 
    mutate(
      fragment_start = round(fragment_start, 0),
      fragment_end   = round(fragment_end, 0),
      fragment_id = paste(fragment_chr, fragment_start, fragment_end, allele, peak_id, sep=":")
    )
  fasta_cell <- file(paste0("fasta_test_",cell_id,".txt"), open = "w")
  for (f in seq_along(1:nrow(cell_frags))){
    chrom = cell_frags[f,"chr"]
    start = cell_frags[f,"fragment_start"] 
    l = cell_frags[f,"fragment_end"]-cell_frags[f,"fragment_start"] 
    fragment_id = cell_frags[f,"fragment_id"]
    fragments_available_alleles <-  cell_frags[f,"available_alleles"]
    for (allele in fragments_available_alleles[[1]]){
      fragment = label_tour$value$genome$get_fragment(chrom, allele,start , l)
      fragment_seq = fragment$get_sequence()
      fragment_id = paste0(">",cell_id,":",cell_frags[f,"fragment_id"],":allele_",allele)
      message(fragment_id,"\n",fragment_seq)
      writeLines(fragment_id, fasta_cell)
      writeLines(fragment_seq, fasta_cell)
    }
  }
  close(fasta_cell)
  i=i+1
  print(i)
  label_tour$step()
}

