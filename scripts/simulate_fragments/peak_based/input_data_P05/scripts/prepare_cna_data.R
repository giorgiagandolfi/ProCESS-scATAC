library(dplyr)
library(data.table)
library(GenomicRanges)
library(readxl)
library(ggplot2)
pat05_cna_cleaned =read.table("../data/DICE_HIPEC_P05_TumourCells_noATACQC_CNA_ConsistentSegmentation_Tidy_Clades_18-08-25.tsv",header = T,sep='\t')
seg_classification =pat05_cna_cleaned %>% 
  filter(CELL!="Diploid") %>%
  filter(CELL!="Clade5") %>%
  filter(CN.states!=2) %>% 
  mutate(seg_id=paste(chrom,start,end,sep=':')) %>% 
  group_by(seg_id) %>% 
  mutate(prop_clade=n()/4) %>% 
  mutate(seg_type=case_when(prop_clade==1~'clonal',
                            TRUE~'subclonal'))

clonal_seg_merged =seg_classification %>% 
  filter(seg_type=='clonal') %>% 
  filter(CN.states!=2) %>% 
  select(chrom,start,end,prop_clade,CN.states) %>% 
  distinct() %>% 
  arrange(chrom, start) %>%
  group_by(chrom) %>%
  mutate(
    grp = cumsum(
      row_number() == 1 |
        start != lag(end) + 1
    )
  ) %>%
  group_by(chrom, grp) %>%
  summarise(
    start = first(start),
    end = last(end),
    prop_clade = mean(prop_clade),
    CN.states=mean(CN.states),
    .groups = "drop"
  ) %>%
  mutate(seg_id = paste(chrom, start, end, sep = ":")) %>%
  select(seg_id, chrom, start, end, prop_clade,CN.states) %>% 
  filter(!chrom%in%c("chrX","chrY")) %>% 
  mutate(cn_len=end-start) 

###### only the clonal segments
clonal_seg_merged_selected = clonal_seg_merged %>% 
  filter(cn_len>=1e7) %>% 
  filter(chrom%in%c("chr5",'chr11','chr17')) %>% 
  mutate(CN.states=round(CN.states,0)) %>% 
  mutate(type=case_when(CN.states==1~"D",
                        TRUE~"A")) %>% 
  mutate(genetic_clone='0')

subclonal_seg_merged =seg_classification %>% 
  filter(seg_type!='clonal') %>% 
  filter(CN.states!=2) %>% 
  select(chrom,start,end,prop_clade,CN.states) %>% 
  distinct() %>% 
  arrange(chrom, start) %>%
  group_by(chrom) %>%
  mutate(
    grp = cumsum(
      row_number() == 1 |
        start != lag(end) + 1
    )
  ) %>%
  group_by(chrom, grp) %>%
  summarise(
    start = first(start),
    end = last(end),
    prop_clade = mean(prop_clade),
    CN.states=mean(CN.states),
    .groups = "drop"
  ) %>%
  mutate(seg_id = paste(chrom, start, end, sep = ":")) %>%
  select(seg_id, chrom, start, end, prop_clade,CN.states) %>% 
  filter(!chrom%in%c("chrX","chrY"))

###### CLONE 1 private CN events
private_cn_clone_1 = seg_classification %>% 
  filter(seg_type!='clonal') %>% 
  filter(CN.states!=2) %>% 
  filter(CELL=="Clade1") %>% 
  filter(prop_clade==0.25) %>% 
  arrange(chrom, start) %>%
  group_by(chrom) %>%
  mutate(
    grp = cumsum(
      row_number() == 1 |
        start != lag(end) + 1
    )
  ) %>%
  group_by(chrom, grp) %>%
  summarise(
    start = first(start),
    end = last(end),
    prop_clade = mean(prop_clade),
    CN.states=mean(CN.states),
    .groups = "drop"
  ) %>%
  mutate(seg_id = paste(chrom, start, end, sep = ":"))

private_cn_clone_1_selected = private_cn_clone_1%>% 
  mutate(cn_len=end-start) %>% 
  filter(cn_len>=1e7) %>%
  select(!grp) %>% 
  distinct() %>% 
  mutate(CN.states=round(CN.states,0)) %>% 
  mutate(type=case_when(CN.states==1~"D",
                        TRUE~"A")) %>% 
  mutate(genetic_clone='A')



###### CLONE 2 private CN events
# use the MYC gene as amplified

###### CLONE 3-4 private CN events
private_cn_clone_34 = seg_classification %>% 
  filter(seg_type!='clonal') %>% 
  filter(CN.states!=2) %>% 
  filter(CELL%in%c("Clade3","Clade4")) %>% 
  filter(prop_clade==0.50) %>% 
  arrange(chrom, start) %>%
  group_by(chrom) %>%
  mutate(
    grp = cumsum(
      row_number() == 1 |
        start != lag(end) + 1
    )
  ) %>%
  group_by(chrom, grp) %>%
  summarise(
    start = first(start),
    end = last(end),
    prop_clade = mean(prop_clade),
    CN.states=mean(CN.states),
    .groups = "drop"
  ) %>%
  mutate(seg_id = paste(chrom, start, end, sep = ":"))

private_cn_clone_34_selected = private_cn_clone_34%>% 
  mutate(cn_len=end-start) %>% 
  filter(cn_len>=1e7) %>%
  mutate(CN.states=round(CN.states,0)) %>% 
  select(!grp) %>% 
  distinct() %>% 
  mutate(type=case_when(CN.states==1~"D",
                        TRUE~"A")) %>% 
  mutate(genetic_clone='C')



###### CLONE 4 private CN events
private_cn_clone_4 = seg_classification %>% 
  filter(seg_type!='clonal') %>% 
  filter(CN.states!=2) %>% 
  filter(CELL%in%c("Clade4")) %>% 
  filter(prop_clade==0.25) %>% 
  arrange(chrom, start) %>%
  group_by(chrom) %>%
  mutate(
    grp = cumsum(
      row_number() == 1 |
        start != lag(end) + 1
    )
  ) %>%
  group_by(chrom, grp) %>%
  summarise(
    start = first(start),
    end = last(end),
    prop_clade = mean(prop_clade),
    CN.states=mean(CN.states),
    .groups = "drop"
  ) %>%
  mutate(seg_id = paste(chrom, start, end, sep = ":"))

private_cn_clone_4_selected = private_cn_clone_4%>% 
  mutate(cn_len=end-start) %>% 
  filter(cn_len>=1e7) %>%
  mutate(CN.states=round(CN.states,0)) %>% 
  select(!grp) %>% 
  distinct() %>% 
  mutate(type=case_when(CN.states==1~"D",
                        TRUE~"A")) %>% 
  mutate(genetic_clone='D')


final_pat05_cn <- do.call("rbind",list(clonal_seg_merged_selected,private_cn_clone_1_selected,private_cn_clone_34_selected,private_cn_clone_4_selected))
saveRDS(object = final_pat05_cn,"../data/final_pat05_copy_number.rds")

### epigenetic classification of cells
epi_cells = read.csv(file = "../data/magic_gsva_bootstrap_zscoreParam_consensus_clusters.csv") %>% 
  tidyr::separate(cell_id,into = c("sample_id","cell_id"),sep = '#',remove = T)
### genetic classification of cells
genetic_cells = read.csv(file = "../data/RandomForest_PredictedClades.csv")

epig_genetic_clades = full_join(epi_cells,genetic_cells)
epig_genetic_clades %>% 
  filter(!is.na(Clade_pred)) %>% 
  filter(!is.na(consensus_cluster)) %>% 
  mutate(consensus_cluster=as.factor(consensus_cluster)) %>% 
  group_by(Clade_pred,consensus_cluster) %>% 
  summarise(n=n()) %>% 
  ggplot(aes(x=Clade_pred,y=n,fill=consensus_cluster,group=consensus_cluster))+
  # geom_bar(stat="identity")
  geom_bar(position="fill", stat="identity")

