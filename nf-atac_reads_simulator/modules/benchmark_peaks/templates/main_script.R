#!/usr/bin/env Rscript
# library(ProCESS)
library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(ggplot2)
library(patchwork)

source("${projectDir}/bin/utils.R")
called_peaks <- read.table("${peaks}", header = FALSE, sep = "\\t")
colnames(called_peaks) <- c(
  "chrom", "start", "end", "name", "score",
  "strand", "signalValue", "pValue", "qValue", "peak"
)
# All per-cell RDS files (grouped by sample_id, staged as a list)
#cell_files <- list.files(path = ".", pattern = "fragments_cell_id_.*", full.names = TRUE)
x <- "${cell_rds_files}"
files <- strsplit(x, " ")[[1]]
peak_accessibility_list <- lapply(files, readRDS) %>% bind_rows()

compare_res = evaluate_peaks(truth = truth_peaks,called = called_peaks,min_overlap = 10)
compare_res = evaluate_peaks_new(truth = truth_peaks,called = called_peaks,min_overlap = 10)

all_peaks_classification <- compare_res\$peak_classification_df
all_peaks_classification <- all_peaks_classification %>% 
  dplyr::mutate(peak_bulk_cell_class=case_when(simulated_pct_cells>=0.8~'high bulk',
                                        simulated_pct_cells<0.8 & simulated_pct_cells>+0.4~'medium bulk',
                                        TRUE~'low bulk'))

all_peaks_classification %>% 
  dplyr::group_by(peak_bulk_cell_class,status) %>% 
  dplyr::summarise(n=dplyr::n())

matching_df=compare_res\$peak_matching_df
matching_df <-matching_df %>% 
  dplyr::mutate(base_difference_start=(simulated_start-inferred_start)) %>% 
  dplyr::mutate(base_difference_end=(simulated_end-inferred_end)) %>% 
  dplyr::mutate(abs_diff=abs(base_difference_end)+abs(base_difference_start))

abs_diff_plt=matching_df %>% 
  ggplot(aes(y=abs_diff,x=""))+geom_boxplot(outliers = F)+
  theme_minimal()
