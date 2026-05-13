set.seed(0)
library(ggplot2)
library(purrr)
library(tidyverse)
library(Seurat)
library(Signac)
library(Rsamtools)
library(GenomicRanges)
library(MASS)

dir_ATAC <- file.path("/data/rds/DMP/UCEC/GENEVOD/ggandolfi/HTAN_data/snATACseq/level_3/")
setwd(dir_ATAC)
## get data from the primary tumour

# barcodes_primary <- readLines("CM618C1-S1-barcodes.tsv")
# peaks_primary <- read.table("CM618C1-S1-peaks.bed", sep="\t")
# peaks_primary_gr <- GRanges(
#   seqnames = peaks_primary$V1,
#   ranges = IRanges(start = peaks_primary$V2, end = peaks_primary$V3)
# )
# peaknames_primary <- paste(peaks_primary$V1, peaks_primary$V2, peaks_primary$V3, sep="-")
# fragpath_primary <- 'CM618C1-S1-fragments.tsv.gz'
process_primary <- readRDS("../level_4/CRC_CM618C2-S1Y2.rds")
metadata = process_primary@meta.data %>% 
  dplyr::rename(barcode=Original_barcode)
barcodes_process_primary <- process_primary$Original_barcode
frags_prim <- CreateFragmentObject(path = fragpath_primary, cells = barcodes_process_primary)

fragments_chr22=read.table(file = 'CM618C1-S1-fragments-top1000.tsv',header = F,sep = '\t')
colnames(fragments_chr22) =c('chr','start','end','barcode','readCount')
### select only a specifi cell type
fragments_chr22 = fragments_chr22 %>% 
  left_join(y = metadata,by = 'barcode') %>% 
  filter(!is.na(cell_type))
selected_frags <-fragments_chr22 %>% filter(cell_type=='Tumor')




fragments_chr22=read.table(file = 'CM618C1-S1-fragments-chr22.tsv',header = F,sep = '\t')
colnames(fragments_chr22) =c('chr','start','end','barcode','readCount')
### select only a specifi cell type
fragments_chr22 = fragments_chr22 %>% 
  left_join(y = metadata,by = 'barcode') %>% 
  filter(!is.na(cell_type))
selected_frags <-fragments_chr22 %>% filter(cell_type=='Tumor')
selected_frags %>% mutate(fragm_len=end-start) %>% ggplot(aes(x=fragm_len))+geom_histogram(binwidth = 1)
# all_pos =selected_frags %>% 
#   # filter(cell_type%in%c("Endothelial",'Tumor')) %>% 
#   ggplot(aes(x=start,fill=cell_type))+
#   geom_histogram()+
#   theme_minimal()
# 
# 
# 
# 
# 
# selected_frags = selected_frags %>% 
#   mutate(fr_len=end-start)
# peak_1 = peaks_primary %>% head(1)
# colnames(peak_1) = c('chr','start','end')
# cuts_real_data = c(selected_frags$start,selected_frags$end)
# 
# region_size = 1e6
# region_start = 1e8
# region_end = region_start+region_size
# 
# selected_frags %>% 
#   filter(start>=region_start) %>% 
#   filter(end<=region_end) %>% 
#   group_by(barcode) %>% 
#   summarise(n_tn5=n())
# 
# 
# # your dataframe
# dt <- as.data.table(selected_frags)
# 
# # -----------------------------------
# # 1. Create fragment midpoint
# # -----------------------------------
# dt[, midpoint := floor((start + end) / 2)]
# 
# # -----------------------------------
# # 2. Define bin size
# # -----------------------------------
# bin_size <- 1000   # 1 kb bins
# 
# # -----------------------------------
# # 3. Assign bins
# # -----------------------------------
# dt[, bin_start := floor(midpoint / bin_size) * bin_size]
# dt[, bin_end := bin_start + bin_size]
# 
# # -----------------------------------
# # 4. Count fragments per bin
# # -----------------------------------
# bin_counts <- dt[
#   ,
#   .(
#     fragments = .N,
#     total_reads = sum(readCount)
#   ),
#   by = .(chr, bin_start, bin_end)
# ]
# 
# # sort by highest signal
# bin_counts <- bin_counts[order(-fragments)]
# 
# 
# # -----------------------------------
# # 1. Create fragment midpoint
# # -----------------------------------
# dt[, midpoint := floor((start + end) / 2)]
# 
# # -----------------------------------
# # 2. Define bin size
# # -----------------------------------
# bin_size <- 1000   # 1 kb bins
# 
# # -----------------------------------
# # 3. Assign bins
# # -----------------------------------
# dt[, bin_start := floor(midpoint / bin_size) * bin_size]
# dt[, bin_end := bin_start + bin_size]
# 
# # -----------------------------------
# # 4. Count fragments per bin
# # -----------------------------------
# bin_counts <- dt[
#   ,
#   .(
#     fragments = .N,
#     total_reads = sum(readCount)
#   ),
#   by = .(chr, bin_start, bin_end)
# ]
# 
# # sort by highest signal
# bin_counts <- bin_counts[order(-fragments)]
# 
# 
# 
# # cuts_real_data = selected_frags$start
# frag_len_dist = selected_frags$fr_len
# 
# #### fit distirbution to fragment lenghts
# fit <- fitdistr(frag_len_dist, densfun = "poisson") 
# df_fragm = tibble(
#   id = seq_along(1:length(cuts_real_data)),
#   cut_site=cuts_real_data
# )
# # 
# # 
# # 
# # hist(frag_len_dist, prob = TRUE, main = "Normal Fit")
# # curve(dnorm(x, mean = fit$estimate[1], sd = fit$estimate[1]), 
# #       add = TRUE, col = "red", lwd = 2)
# 
# 
# 
# # median number of cuts per cell in a region
# selected_frags %>% 
#   ggplot(aes(x=fr_len))+
#   geom_histogram(binwidth = 1)
# 
# 
# 
# # numero di cellule nel campione
# num_of_cells <- 2000
# # lunghezza della regione investigata (regioni aperte + chiuse)
# region_size <- max(selected_frags$end)-min(selected_frags$start)
# 
# 
# 
# df_sampled_cuts = tibble(
#   id = NA,#seq_along(1:num_of_cells),
#   fr_start = 0,
#   fr_end=0,
#   fr_len=0
# )
# 
# for (i in seq_along(1:num_of_cells)) {
#   start_cut =sample(x = cuts_real_data,size = 1)
#   fragm_len = sample(x = frag_len_dist,size = 1)
#   end_cut=start_cut+fragm_len
#   df_sampled_cuts = rbind(df_sampled_cuts,tibble(
#     id = paste0('cell_',i),
#     fr_start = start_cut,
#     fr_end=end_cut,
#     fr_len=fragm_len
#   ))
# }
# df_sampled_cuts = df_sampled_cuts %>% 
#   filter(!is.na(id))
# df_sampled_cuts %>% 
#   ggplot(aes(x=fr_len))+
#   geom_histogram(binwidth = 1)
# 
#                