library(ggplot2)
library(dplyr)
library(extraDistr)
library(readxl)
library(tidyverse)
library(Seurat)
library(Signac)
library(GenomicRanges)
dir_ATAC <- file.path("/data/rds/DMP/UCEC/GENEVOD/ggandolfi/HTAN_data/snATACseq/")

source('/data/rds/DMP/UCEC/GENEVOD/ggandolfi/Github/ProCESS-scATAC/scripts/utils_genomics.R')
setwd(dir_ATAC)
genome_coordinates = readRDS('level_3/genome_coordinates_hg38.rds')
chr1_coordinates =genome_coordinates %>% filter(chr=='chr1') %>% data.table()

bin_size <- 1e4   # 1 M

binned_chr1 =binRegion(start = chr1_coordinates$from,end = chr1_coordinates$to,binSize = bin_size) %>% as.data.frame()
binned_chr1 <- binned_chr22 %>%
  mutate(
    abs_start = chr1_coordinates$from + start,
    abs_end   = chr1_coordinates$from + end
  )
fragments_chr22=read.table(file = file.path(dir_ATAC,'level_3','CM618C1-S1-fragments-top1000.tsv'),header = F,sep = '\t')
colnames(fragments_chr22) =c('chr','start','end','barcode','readCount')
### select only a specifi cell type
metadata=tumor_subset@meta.data %>% 
  dplyr::rename(barcode=Original_barcode)
fragments_chr22 = fragments_chr22 %>% 
  left_join(y = metadata,by = 'barcode') %>% 
  filter(!is.na(cell_type))
selected_frags <-fragments_chr22 %>% filter(cell_type=='Tumor') %>% 
  mutate(from = start + vfrom[chr],
         to = end + vfrom[chr]) %>% 
  mutate(fragment_len=to-from) 

cut_sites <- selected_frags %>%
  dplyr::select(chr, from, to, barcode) %>%
  pivot_longer(
    cols = c(from, to),
    names_to = "type",
    values_to = "cut_site"
  ) %>%
  dplyr::select(chr, cut_site, barcode)




# cutting_sites = c(selected_frags$start,selected_frags$end)
# cutting_sites = sort(cutting_sites)


cuts_in_bins <- cut_sites %>%
  mutate(
    binID = findInterval(cut_site, binned_chr1$abs_start)
  ) %>%
  filter(
    cut_site <= binned_chr1$abs_end[binID]
  )

cuts_in_bins %>% 
  ggplot(aes(x=cut_site))+geom_histogram()+
  theme_bw()+
  xlab("chr1")+
  ylab('N of cutting sites')

# count cuts per barcode per bin
barcode_counts <- cuts_in_bins %>%
  group_by(binID, barcode) %>%
  summarise(n_cuts = n(), .groups = "drop")


bin_medians <- barcode_counts %>%
  group_by(binID) %>%
  summarise(
    median_cuts = median(n_cuts),
    mean_cuts = mean(n_cuts),
    sd_cuts = sd(n_cuts),
    .groups = "drop"
  ) %>%
  full_join(binned_chr1, by = "binID")
bin_size <- 200

high_prob_regions <- cut_sites %>%
  mutate(
    bin = floor(cut_site / bin_size) * bin_size
  ) %>%
  group_by(chr, bin) %>%
  mutate(
    n_cuts = n()
    # .groups = "drop"
  ) %>%
  arrange(desc(n_cuts))




high_prob_regions <- high_prob_regions %>%
  mutate(
    start = bin,
    end = bin + bin_size - 1
  ) %>% 
  filter(n_cuts>=20)

library(IRanges)

# build ranges for bins
bins_ir <- IRanges(
  start = bin_medians$abs_start,
  end   = bin_medians$abs_end
)

# build ranges for high prob regions
regions_ir <- IRanges(
  start = high_prob_regions$start,
  end   = high_prob_regions$end
)

# find overlaps
hits <- findOverlaps(regions_ir, bins_ir)

mapped <- tibble(
  region_id = queryHits(hits),
  bin_id    = subjectHits(hits)
)

# attach metadata
result <- mapped %>%
  left_join(high_prob_regions %>% mutate(region_id = row_number()), by = "region_id") %>%
  left_join(bin_medians %>% mutate(bin_id = row_number()), by = "bin_id", suffix=c('_acc_reg','_genome_bin'))

result <- result %>% 
  dplyr::select(chr,cut_site,barcode,start_acc_reg,end_acc_reg,median_cuts,sd_cuts,start_genome_bin,end_genome_bin)

#simualte number of cells
num_of_cells = 1000
selected_end_genome_bin = result[1,('end_genome_bin')]
selected_start_genome_bin = result[1,('start_genome_bin')]
region_size= 1e4
open_regions =result %>% 
  dplyr::filter(start_genome_bin==selected_start_genome_bin) %>% 
  dplyr::filter(end_genome_bin==selected_end_genome_bin) %>% 
  dplyr::select(start_acc_reg,end_acc_reg,start_genome_bin,end_genome_bin) %>%
  dplyr::filter(!is.na(start_acc_reg)) %>% 
  unique() %>% 
  ungroup()

open_regions <- open_regions %>%
  mutate(
    new_start_acc = 1 +
      (start_acc_reg - start_genome_bin) *
      (region_size - 1) /
      (end_genome_bin - start_genome_bin),
    
    new_end_acc = 1 +
      (end_acc_reg - start_genome_bin) *
      (region_size - 1) /
      (end_genome_bin - start_genome_bin)
  )
# open_regions = open_regions %>% head(1)
univocity_read_size <- 10
univocity_read_stddev <- 5


C <- array(0, dim = c(region_size))
mask <- array(FALSE, dim = c(region_size))

for (region in 1:nrow(open_regions)) {
  print(region)
  start_reg = as.numeric(open_regions[region,'new_start_acc'])
  end_reg = as.numeric(open_regions[region,'new_end_acc'])
  for (i in start_reg:end_reg) {
    mask[i] <- TRUE
  }
}
df <- data.frame(
  index = 1:length(C),
  value = C,
  cell=NA,
  open = mask
)
cut_tn5 <- data.frame(
  index = NA,
  cut_site = NA,
  cell =NA
)
fragments <- tibble(
  start = 0,
  end   = 0,
  cell =NA
)

for (i in seq(1, num_of_cells)) {
  df_tmp <- data.frame(
    index = 1:length(C),
    value = C,
    cell=paste0("cell_",i),
    open = mask
  )
  num_of_cuts <- 4
  # num_of_cuts <- as.integer(rnorm(1, mean = num_of_cuts_mean,
  #                                 sd = num_of_cuts_stdev))
  
  cuts <- c()
  for (j in 1:num_of_cuts) {
    
    cut_set <- FALSE
    
    while (!cut_set) {
      cut_pos <- sample(1:length(C), 1)
      cut_set <- mask[cut_pos]
    }
    
    cuts <- c(cuts, cut_pos)
  }
  
  cuts <- sort(cuts)
  cut_tn5 <- rbind(cut_tn5,data.frame(
    index = 1:length(cuts),
    cut_site = cuts,
    cell = paste0("cell_",i)
  ))
  j <- 1
  while (j <= length(cuts)) {
    if (j == length(cuts)) {
      #cut_size <- size-cuts[j]-1
    } else {
      cut_size <- cuts[j+1]-cuts[j]
    }
    if ((cut_size > 2*univocity_read_size |
        cut_size>=rnorm(1, mean = 2*univocity_read_size, sd = 5)) & cut_size<=400) {
      C[cut_size] <- C[cut_size]+1
      fragments <- rbind(fragments,tibble(
        start = cuts[j],
        end   = cuts[j+1],
        cell = paste0("cell_",i)
      ))
    }
    j <- j+1
  }
  fragments <- na.omit(fragments)
  df <- rbind(df,df_tmp)
  print(i)
}

df <- df %>% 
  filter(!is.na(cell))

cut_tn5 <- cut_tn5 %>% filter(!is.na(index))
cut_tn5 <- cut_tn5 %>%
  mutate(cell_y = as.numeric(factor(cell)))
fragments <- fragments %>%
  mutate(cell_y = as.numeric(factor(cell)))


open_df <- open_regions %>%
  mutate(status = "open") %>% 
  dplyr::select(new_start_acc,new_end_acc,status) %>% 
  dplyr::rename(end=new_end_acc) %>% 
  dplyr::rename(start=new_start_acc)

# closed regions (gaps)
closed_df <- tibble(
  start = c(1, head(open_df$end, -1) + 1),
  end   = c(head(open_df$start, -1) - 1,
            region_size)
) %>%
  filter(start <= end) %>%
  mutate(status = "closed")

# final dataframe
chromatin_acc_df <- bind_rows(open_df, closed_df) %>%
  arrange(start)

fragments <- fragments %>% 
  mutate(fragm_len=end-start) %>% 
  mutate(fragm_class=case_when(fragm_len<147 ~ "NF",
                               fragm_len>=147 & fragm_len<250 ~ "mono",
                               TRUE ~ "bi/tri" ))
ggplot() +
  
  # chromatin accessibility annotation
  geom_segment(
    data = chromatin_acc_df, aes(x = start,xend = end,y = 0,yend = 0, color = status),
    linewidth = 3
  ) +
  
  # fragments
  geom_segment(
    # data =fragments_sampled,
    data = fragments %>% filter(start != 0 & end != 0) %>% head(20),
    aes(x = start,xend = end,y = cell_y,yend = cell_y),
    linewidth = 2,alpha = 0.5) +
  # cut sites
  # geom_point(data = cut_tn5,aes(x = cut_site,y = cell_y),
  #   size = 0.2,
  #   alpha = 0.5
  # ) +
  theme_minimal() +
  # facet_wrap(~fragm_class,scales = "free")+
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )


fragments %>% 
  ggplot(aes(x=fragm_len))+
  geom_histogram(binwidth = 1)+
  theme_bw()

