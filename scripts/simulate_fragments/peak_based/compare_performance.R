library(dplyr)
library(ggplot2)
library(tidyverse)
library(ComplexHeatmap)
# library(circlize)
library(data.table)
library(readxl)

selected_frags_dist_len = readRDS("selected_frags_dist_len.rds")
df_peak_final = readRDS("df_peak_final.rds")
df_peak_final <- df_peak_final %>% 
  separate(col = peak_id,into = c("chr","from","to"),sep = ":",remove = F) %>% 
  mutate(from=as.numeric(from),
         to=as.numeric(to))
genome_coordinates = readRDS('genome_coordinates_hg38.rds')


chr_start <- genome_coordinates$from[genome_coordinates$chr == "chr22"]
df_peak_final$start_rel <- df_peak_final$from - chr_start + 1
df_peak_final$end_rel   <- df_peak_final$to   - chr_start + 1
df_peak_final = df_peak_final %>% dplyr::select(chr,start_rel,end_rel,peak_id,cell_id,status)

macs_peaks <- read.table("test_peaks.narrowPeak", header = F, sep = "\t")
colnames(macs_peaks) <- c(
  "chrom", "start", "end", "name", "score",
  "strand", "signalValue", "pValue", "qValue", "peak"
)

library(GenomicRanges)

df_peak_final_cell <- df_peak_final %>% 
  filter(cell_id==92086) %>% 
  filter(status==1)
gt <- GRanges(
  seqnames = gsub(pattern = "chr",replacement = "",x = df_peak_final_cell$chr),
  ranges = IRanges(start = df_peak_final_cell$start_rel, end = df_peak_final_cell$end_rel))

called <- GRanges(
  seqnames = macs_peaks$chrom,
  ranges = IRanges(start = macs_peaks$start, end = macs_peaks$end)
)

hits <- findOverlaps(gt,called)
queryHits(hits)
subjectHits(hits)
matching_peaks = cbind(df_peak_final_cell[queryHits(hits),],macs_peaks[subjectHits(hits),])
matching_peaks = matching_peaks %>% 
  dplyr::rename(chr_gt=chr,
                start_gt=start_rel,
                end_gt=end_rel,
                chr_macs=chrom,
                start_macs=start,
                end_macs=end
                )
matching_peaks = matching_peaks %>% 
  mutate(abs_diff_start=abs(start_gt-start_macs)) %>% 
  mutate(abs_diff_end=abs(end_gt-end_macs)) %>% 
  mutate(overall_offset=abs_diff_end+abs_diff_start)


overallp_performance = full_join(df_peak_final_cell,matching_peaks,by="peak_id")
df_peak_final_cell_filtered = df_peak_final_cell[queryHits(hits),c("chr","start_rel","end_rel")]
df_peak_final_cell_filtered = df_peak_final_cell_filtered %>% 
  mutate(type="ground truth") %>% 
  dplyr::rename(chrom=chr,
                start=start_rel,
                end=end_rel)  %>% 
  mutate(peak_id=paste0("peak_id",seq_along(1:nrow(df_peak_final_cell_filtered))))

macs_peaks_filtered = macs_peaks[subjectHits(hits),c("chrom","start","end")]
macs_peaks_filtered = macs_peaks_filtered %>% 
  mutate(type="inferred") %>% 
  mutate(peak_id=paste0("peak_id",seq_along(1:nrow(macs_peaks_filtered))))

compare_peaks = rbind(df_peak_final_cell_filtered,macs_peaks_filtered)
compare_peaks %>% 
  # filter(peak_id=="peak_id1") %>%
  ggplot(aes(x=start,xend=end,y=peak_id,yend=peak_id,color=type)) +
  geom_segment()+
  facet_wrap(~peak_id,scales = "free",ncol = 1)+
  theme(strip.text.x = element_blank(),axis.text.x = element_blank())


