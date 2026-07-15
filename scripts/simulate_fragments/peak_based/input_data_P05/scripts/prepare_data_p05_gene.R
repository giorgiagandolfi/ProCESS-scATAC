library(dplyr)
library(purrr)
library(tidyverse)
library(readxl)
library(GenomicRanges)
source("/data/rds/DMP/UCEC/GENEVOD/ggandolfi/Github/ProCESS-scATAC/scripts/simulate_fragments/peak_based/utils.R")

#######
# This script prepares the data in order to be used later on in ProCESS simulation
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
clonal_peaks = readRDS('/data/rds/DMP/UCEC/GENEVOD/ggandolfi/Github/ProCESS-scATAC/scripts/simulate_fragments/peak_based/input_data_P05/data/a_scores_per_gene_log1p_score.rds')

clonal_peaks = clonal_peaks %>% 
  mutate(peak=paste(seqnames,start,end,sep='-'))

gr_clonal_peaks = GRanges(
  seqnames = clonal_peaks$seqnames,
  ranges = IRanges(start = clonal_peaks$start, end = clonal_peaks$end),
  score = clonal_peaks$log1p_score
)

gr_tissue_peaks = GRanges(
  seqnames = crc_peaks_top_fch$seqnames,
  ranges = IRanges(start = crc_peaks_top_fch$start, end = crc_peaks_top_fch$end)
)


overlaps_tissue_clonal <- findOverlaps(query = gr_clonal_peaks,subject = gr_tissue_peaks,minoverlap = 100)
# overlaps_tissue_changing <- findOverlaps(query = gr_chaning_peaks,subject = gr_tissue_peaks,minoverlap = 100)
overlaped_peaks_clonal <-  clonal_peaks[(queryHits(overlaps_tissue_clonal)),]
to_remove_tissue <- c(subjectHits(overlaps_tissue_clonal)) %>% unique()
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


epigenetic_peaks <- clonal_peaks %>% 
  dplyr::rename(score=log1p_score,
                class=consensus_cluster) %>% 
  select(peak,seqnames,start,end,class,score,gene) %>%
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



p05_peaks <- epigenetic_peaks
sampled_genes = p05_peaks %>% 
  pull(gene) %>% unique() %>%  sample(size = 100)
p05_peaks_sampled = p05_peaks %>% 
  filter(gene%in%c(sampled_genes,"MYC"))

activity_df <- p05_peaks_sampled %>% 
  dplyr::rename(a_score=score) %>% 
  dplyr::select(a_score,gene,mutant,epistate) %>% 
  distinct()

#### create the activity list expected by ProCESS
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
    names(pathways_vect) <- pathways_state$gene
    mutant_activity_list[[epi]]<-pathways_vect
  }
  activity_list[[mut]]<-mutant_activity_list
}
saveRDS(object = activity_list,file = '../data/activity_list_gene_level_log1pscore.rds')



# p05_peaks <- p05_peaks %>% 
#   select(peak,gene) %>% 
#   distinct()


p05_peaks_sampled <- p05_peaks_sampled %>% 
  select(peak,gene) %>% 
  distinct()


saveRDS(object = p05_peaks_sampled,file = '../data/all_peaks_gene_df.rds')

