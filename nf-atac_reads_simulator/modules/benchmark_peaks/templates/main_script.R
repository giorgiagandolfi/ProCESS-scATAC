#!/usr/bin/env Rscript
library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(ggplot2)

source("${projectDir}/bin/utils.R")
called_peaks <- read.table("${peaks}", header = FALSE, sep = "\\t")
colnames(called_peaks) <- c(
  "chrom", "start", "end", "name", "score",
  "strand", "signalValue", "pValue", "qValue", "peak"
)

frags_peak_mapping_files <- strsplit("${cell_peak_txt_files}", " ")[[1]]
simulated_frags_files <- strsplit("${cell_frags_rds_files}", " ")[[1]]

sequenced_peaks_list<-list()
for (i in 1:length(simulated_frags_files)){
  peak_frags_pre_noise = read.table(frags_peak_mapping_files[i])
  peak_frags_pre_noise = peak_frags_pre_noise %>% 
    dplyr::mutate(fragment_start=round(fragment_start,0)) %>% 
    dplyr::mutate(fragment_end=round(fragment_end,0))
  sequenced_frags = readRDS(simulated_frags_files[i])
  sequenced_frags = sequenced_frags %>% 
    dplyr::mutate(fragment_start=round(fragment_start,0)) %>% 
    dplyr::mutate(fragment_end=round(fragment_end,0))
  sequenced_frags = sequenced_frags %>% 
    filter(sequenced==1) %>% 
    filter(fragment_type=='peak') %>% 
    inner_join(peak_frags_pre_noise,by=c('fragment_start','fragment_end','fragment_allele','fragment_chr'))
  
  sequenced_peaks_list[[i]] = sequenced_frags %>% 
    select(peak, peak_chr,peak_from,peak_to,cell_id)
  print(i)
}

sequenced_peaks_df <- do.call('rbind',sequenced_peaks_list)
tot_cells = length(frags_peak_mapping_files)
sequenced_peaks_df=sequenced_peaks_df %>% 
  group_by(peak) %>%
  dplyr::mutate(pct_cells=n()/tot_cells) %>% 
  dplyr::select(!cell_id) %>% 
  dplyr::distinct()


truth=sequenced_peaks_df
called=called_peaks
min_overlap=10
results = evaluate_peaks_new(truth,called,min_overlap)
saveRDS(object = results,file = paste0("${meta.sample_id}","_metrics.rds"))
