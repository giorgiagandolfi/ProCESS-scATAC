# rm(list=ls())
library(dplyr)
library(purrr)
library(tidyverse)
library(readxl)
library(GenomicRanges)
source("utils.R")



####### 

crc_peaks <- read_excel('/data/rds/DMP/UCEC/GENEVOD/ggandolfi/HTAN_data/41586_2023_6682_MOESM4_ESM.xlsx',sheet = 1) %>% 
  dplyr::filter(Cancer=="CRC") %>% 
  dplyr::filter(!grepl("chrX", peak)) %>% 
  dplyr::filter(!grepl("chrY", peak))
crc_peaks_top_fch <- read_excel('/data/rds/DMP/UCEC/GENEVOD/ggandolfi/HTAN_data/41586_2023_6682_MOESM4_ESM.xlsx',sheet = 3) %>% 
  filter(Cancer=="CRC") %>% 
  filter(!grepl("chrX", peak)) %>% 
  filter(!grepl("chrY", peak))
crc_peaks_not_fch = anti_join(crc_peaks,crc_peaks_top_fch)
crc_peaks_top_fch = crc_peaks_top_fch %>% 
  separate(col = peak,into = c('seqnames','start','end'),sep = '-',remove = F) %>% 
  mutate(start=as.numeric(start),
         end=as.numeric(end))
clonal_peaks = readRDS('input_data_P05/peak_per_pathways.rds')
clonal_peaks = clonal_peaks %>% 
  mutate(peak=paste(seqnames,start,end,sep='-'))
# changing_peaks=readRDS('input_data_P05/changing_peaks.rds')
# changing_peaks = changing_peaks %>% 
#   mutate(peak=paste(seqnames,start,end,sep='-'))
gr_clonal_peaks = GRanges(
  seqnames = clonal_peaks$seqnames,
  ranges = IRanges(start = clonal_peaks$start, end = clonal_peaks$end),
  score = clonal_peaks$pathway
)

gr_tissue_peaks = GRanges(
  seqnames = crc_peaks_top_fch$seqnames,
  ranges = IRanges(start = crc_peaks_top_fch$start, end = crc_peaks_top_fch$end)
  # score = clonal_peaks$pathway
)

# gr_chaning_peaks = GRanges(
#   seqnames = changing_peaks$seqnames,
#   ranges = IRanges(start = changing_peaks$start, end = changing_peaks$end)
#   # score = clonal_peaks$pathway
# )

overlaps_tissue_clonal <- findOverlaps(query = gr_clonal_peaks,subject = gr_tissue_peaks,minoverlap = 100)
# overlaps_tissue_changing <- findOverlaps(query = gr_chaning_peaks,subject = gr_tissue_peaks,minoverlap = 100)
overlaped_peaks_clonal <-  clonal_peaks[(queryHits(overlaps_tissue_clonal)),]
to_remove_tissue <- c(subjectHits(overlaps_tissue_clonal),subjectHits(overlaps_tissue_changing)) %>% unique()
overlaped_peaks_tissue <-  crc_peaks_top_fch[to_remove_tissue,]

filtered_tissue_peaks <- crc_peaks_top_fch %>% 
  filter(!peak%in%c(overlaped_peaks_tissue$peak))


fixed_peaks <-filtered_tissue_peaks %>% 
  select(peak,seqnames,start,end) %>% 
  mutate(type='fixed',
         mutant=NA,
         epistate=NA,
         score=1) %>% 
  dplyr::rename(chr=seqnames)


pathways_epigenetic_classes = readRDS('input_data_P05/a_scores_scaled_per_group.rds')
pathways_long <- pathways_epigenetic_classes %>%
  as.data.frame() %>%
  rownames_to_column("pathway") %>%
  pivot_longer(
    cols = -pathway,
    names_to = "class",
    values_to = "score"
  )
epigenetic_peaks <- clonal_peaks %>% 
  left_join(pathways_long,relationship = 'many-to-many') %>% 
  select(peak,seqnames,start,end,class,score,pathway) %>%
  dplyr::rename(chr=seqnames) %>% 
  mutate(mutant=case_when(class%in%c(1,3)~"A",
                          TRUE~"B")) %>% 
  mutate(epistate=case_when(class==1~"+",
                            class==3~"-",
                            class==2~"+",
                            class==4~"-",
                          TRUE~NA)) %>% 
  mutate(type='clonal') %>% 
  select(!class)

epigenetic_states <- c("A+","A-","B+","B-")
fixed_peaks_all <- lapply(epigenetic_states, function(epi){
  fixed_peaks %>% 
    mutate(mutant=strsplit(epi, "")[[1]][1],
           epistate=strsplit(epi, "")[[1]][2]) %>% 
    mutate(pathway='CRC_TISSUE')
}) %>% bind_rows()


# fluctuating_peaks <- changing_peaks %>% 
#   select(peak,seqnames,start,end) %>% 
#   dplyr::rename(chr=seqnames) %>% 
#   mutate(type='fluctuating',
#          mutant=NA,
#          epistate=NA,
#          score=0.1)
  
# p05_peaks <- do.call('rbind',list(fixed_peaks,fluctuating_peaks,epigenetic_peaks))
p05_peaks <- do.call('rbind',list(fixed_peaks_all,epigenetic_peaks))

activity_df <- p05_peaks %>% 
  group_by(mutant,epistate,pathway) %>% 
  summarise(a_score=mean(score))

mutants <- activity_df$mutant %>% unique()
epistates <- activity_df$epistate %>% unique()

activity_list <- list()
for (mut in mutants){
  mutant_activity_list <- list()
  for (epi in epistates){
    pathways_state =activity_df %>% 
      filter(mutant==mut,
             epistate==epi) 
    pathways_vect <- pathways_state$a_score
    names(pathways_vect) <- pathways_state$pathway
    mutant_activity_list[[epi]]<-pathways_vect
  }
  activity_list[[mut]]<-mutant_activity_list
}
saveRDS(object = activity_list,file = 'input_data_P05/activity_list.rds')

p05_peaks <- p05_peaks %>% 
  select(peak,pathway) %>% 
  distinct()
saveRDS(object = p05_peaks,file = 'input_data_P05/all_peaks_df.rds')
p05_peaks %>% pull(peak) %>% unique() %>% length()

p05_peaks %>% 
  select(peak,type) %>% 
  distinct() %>% 
  ggplot(aes(x=type,fill=type)) +
  geom_bar()+
  # scale_fill_manual(values=mutant_cols)+
  theme_minimal()+
  ggtitle(label = "Peak distribution per class")
