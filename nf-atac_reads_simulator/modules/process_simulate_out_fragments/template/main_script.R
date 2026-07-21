#!/usr/bin/env Rscript

# Template R script for in-peak fragment simulation
# Nextflow template variables are injected at runtime

# --- Inputs from Nextflow process ---
cell_id<- "${meta.id}"
sample_id<- "${meta.sample_id}"
in_peak_fragments_file<- "${in_peak_frg}"
gap_file <- "${gap_file}"
centromere_file <- "${centromere_file}"
fragm_len_dist_out <- "${frag_dist_out}"
prefix<- "${prefix}"
args<- "${args}"


library(ProCESS)
library(dplyr)
library(tidyverse)
library(parallel)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
source("../../bin/utils.R")


message("Peak fragments generated")
### extract regions that are not peak
cell_bg_regions <- get_background_regions(peak_fragments_df = chromosomes_frags,genome = 'hg38',
                                          gaps_file = gap_file,filter_small_than = 150,
                                          centromeres_file = centromere_file)
final_mapping <- readRDS(fragm_len_dist_out)
frag_len_out_peak = final_mapping %>% 
  filter(region_type=='out peak') %>% 
  # ggplot(aes(x=fragment_len))+geom_density()
  pull(fragment_len)
frag_len_out_peak_dens <- density(frag_len_out_peak,from=100)



background_frg=simulate_background_fragments(background_regions = cell_bg_regions,
                                             lambda_per_kb = 0.1,frag_len_out_peak_dens = frag_len_out_peak_dens)

message("Background fragments generated")
### get cell info into peaks
outidr <- "fragments_cells_big_with_background_01_lambda_sparsity_085_filtered_peaks_tss/"

dir.create(path = outidr)
chromosomes_frags <- chromosomes_frags %>% inner_join(cell_info)
saveRDS(object = chromosomes_frags,file = paste0(outidr,'cell_',selected_cells,'_all_fragments.rds'))
###### create bed file

sampled_cells = cell_info %>% filter(!is.na(sample)) %>% pull(cell_id) %>% unique()


