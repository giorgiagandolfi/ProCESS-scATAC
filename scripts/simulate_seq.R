rm(list=ls())
library(ProCESS)
library(dplyr)
setwd("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC/")
phylo_forest <- load_phylogenetic_forest("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC/phylo_forest_atac3.sff")
chromosomes <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
chromosomes <- c("1","2","12","17","5")
seq_results <- parallel::mclapply(chromosomes, function(c) {
  simulate_seq(phylo_forest, coverage = 100, chromosomes = c, 
               # output_dir = sam_folder, 
               write_SAM = FALSE,
               purity = 1, 
               with_normal_sample = F)
}, mc.cores = parallel::detectCores()) 
# %>% do.call("bind_rows", .)
#seq_results <- simulate_seq(phylo_forest, coverage = 80)
# list.files(sam_folder)

saveRDS(seq_results, "/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC/seq_results3.rds")
print("sequencing ended")
