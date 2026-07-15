library(ggplot2)
library(dplyr)
library(extraDistr)
library(readxl)
library(tidyverse)
library(Seurat)
library(Signac)
library(GenomicRanges)
dir_ATAC <- file.path("/data/rds/DMP/UCEC/GENEVOD/ggandolfi/HTAN_data/snATACseq/")

# peak_dataset <- read_xlsx("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC/41586_2023_6682_MOESM4_ESM.xlsx",sheet = 1)
fragpath_primary <- file.path(dir_ATAC,'level_3','CM618C1-S1-fragments.tsv.gz')
process_primary <- readRDS(file.path(dir_ATAC,'level_4',"CRC_CM618C2-S1Y2.rds"))
barcodes_process_primary <- process_primary$Original_barcode

frags_prim <- CreateFragmentObject(path = fragpath_primary, cells = barcodes_process_primary)
Fragments(process_primary) <- NULL
# attach new fragment object
Fragments(process_primary) <- frags_prim

tumour_cells_barcodes=process_primary@meta.data %>% 
  filter(cell_type=='Tumor') %>% 
  pull(Original_barcode)

tumor_subset <- subset(
  x = process_primary,
  subset = cell_type == "Tumor"
)


peak_data_single_sample <- rownames(tumor_subset@assays$pancan$counts)
peak_data_single_sample <- as.data.frame(peak_data_single_sample)
colnames(peak_data_single_sample) <- "peak"
crc_peaks <- peak_data_single_sample %>% 
  # filter(Cancer=="CRC") %>% 
  separate(col = peak,into = c("chr","from","to"),sep = "-",remove = F) %>% 
  filter(chr%in%c("chr22")) %>% 
  mutate(from=as.numeric(from)) %>% 
  mutate(to=as.numeric(to)) %>% 
  mutate(peak_size=to-from)


genome_coordinates <- readRDS(file.path(dir_ATAC,'level_3','genome_coordinates_hg38.rds'))
vfrom = genome_coordinates$from
names(vfrom) = genome_coordinates$chr

crc_peaks <- crc_peaks %>%
  mutate(from = from + vfrom[chr],
         to = to + vfrom[chr]) %>%
  mutate(peak_id=paste(chr,from,to,sep=":"))
saveRDS(file="crc_peaks.rds",object = crc_peaks)
fragments_chr22=read.table(file = file.path(dir_ATAC,'level_3',"chrom_split","chr8.fragments.tsv"),header = F,sep = '\t')
colnames(fragments_chr22) =c('chr','start','end','barcode','readCount')
selected_frags <-fragments_chr22 %>% 
  mutate(fragment_len=end-start) 
selected_frags_dist_len = selected_frags %>% pull(fragment_len)
ggplot(selected_frags,aes(x=fragment_len))+geom_histogram(binwidth = 1)
saveRDS(object = selected_frags_dist_len,file="selected_frags_dist_len.rds")

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
selected_frags_dist_len = selected_frags %>% pull(fragment_len)
saveRDS(object = selected_frags_dist_len,file="selected_frags_dist_len.rds")
set.seed(123)

n_peaks <- nrow(crc_peaks)
n_peaks <- 100



cna_data_sim = readRDS('/data/rds/DMP/UCEC/GENEVOD/ggandolfi/process_simulation/cell_fragmentation.rds')
#### convert to absolute coordinates
cna_data_sim = cna_data_sim %>% 
  mutate(chr=paste0('chr',chr)) %>% 
  mutate(from = begin + vfrom[chr],
         to = end + vfrom[chr])

n_cells <- cna_data_sim %>% pull(cell_id) %>% unique() %>% length()
n_cells <- 10
cell_ids <- cna_data_sim %>% pull(cell_id) %>% unique() 
# n_peaks <- 20


list_final_frag <- list()
peak_df_with_cna = tibble()
for (c in 1:n_cells){
  print(paste0(c,'/',n_cells))
  list_cell_frag <- list()
  for (p in 1:n_peaks){
    # print(p)
    ### extract peak information
    peak_df = crc_peaks[p,]
    
    
    #### peak cna
    cell_cna = cna_data_sim %>% filter(cell_id==cell_ids[c]) %>% 
      select(chr,major,minor,from,to)
    
    peak_df = peak_df %>%
      left_join(cell_cna, by = "chr") %>%
      filter(
        from.x <= to.y,
        to.x   >= from.y
      ) %>%
      transmute(peak,peak_id,chr,peak_from = from.x,
                peak_to= to.x,major, minor)
    peak_df_with_cna=rbind(peak_df,peak_df_with_cna)
    total_cna_peak=peak_df$major+peak_df$minor
    # print(as.character(crc_peaks[p,"peak_id"]))
    # print(total_cna_peak)
    
    list_cell_frag[[p]] <-sample_fragments_for_peak(peak_id = as.character(crc_peaks[p,"peak_id"]),
                                                     peak_from = as.numeric(crc_peaks[p,"from"]),
                                                    peak_to = as.numeric(crc_peaks[p,"to"]),
                                                    fragment_len_dist = selected_frags_dist_len,
                                                    tot_cn = total_cna_peak) %>% 
      mutate(cellID=cell_ids[c])
  }
  list_final_frag[[c]] <- do.call("rbind",list_cell_frag)
}

results <- do.call("rbind",list_final_frag)


#########################################
results$from <- as.numeric(results$from)
results$to <- as.numeric(results$to)
results<- results %>% 
  mutate(frag_final=paste0(fragment,"_",peak)) %>% 
  mutate(frag_size=to-from)


starting_dist_frag_sizes <- replicate(
  100000,
  sample_fragment_size()
)

library(dplyr)
library(ggplot2)

ggplot() +
  
  # simulated distribution
  geom_histogram(
    data = results,
    aes(x = frag_size, y = after_stat(density)),
    binwidth = 1,
    fill = "steelblue",
    alpha = 0.5
  )+
  
  # input distribution
  # geom_histogram(
  #   data = selected_frags %>% filter(fragment_len <= 800),
  #   aes(x = fragment_len, y = after_stat(density)),
  #   binwidth = 1,
  #   fill = "tomato",
  #   alpha = 0.5
  # ) +
  geom_density(
    data=data.frame(fragment_len=starting_dist_frag_sizes),
    # data = selected_frags %>% filter(fragment_len <= 800),
    aes(x = fragment_len),
    binwidth = 1,
    # fill = "tomato",
    color="tomato",
    alpha = 0.3
  )#+
  # ggtitle(label = paste0('N cells: ',n_cells,'\nN peaks: ',n_peaks))+  theme_minimal()
  

results %>% 
  group_by(peak,cellID) %>% 
  summarise(
    n_fragments = n_distinct(frag_final),
    .groups = "drop"
  ) %>% 
  group_by(n_fragments) %>% 
  summarise(
    n=n(),
    .groups = "drop"
  ) %>% 
  ggplot(aes(x = "", y = n, fill = factor(n_fragments))) +
  geom_col(width = 1,color='white') +
  coord_polar(theta = "y") +
  labs(fill = "Fragments") +
  scale_fill_manual(values=c("1"='royalblue3','2'='steelblue1'))+
  theme_void()

results_wt_cna <- results %>% dplyr::rename(peak_id=peak) %>% left_join(peak_df_with_cna)

results_wt_cna %>% 
  mutate(karyotype=paste0(major,':',minor)) %>% 
  # filter(karyotype!='1:1') %>%
  filter(peak_id=='chr22:2844481288:2844481788') %>%
  # filter(cellID%in%c('5531','2917')) %>% 
  # separate(col = peak,into = c('chr','from_peak','to_peak'),sep = ':',remove = F) %>%
  ggplot()+
  geom_segment(
    aes(
      x = from,
      xend = to,
      y = allele,
      yend = allele, color=allele
    ),
    linewidth = 2
    # color='goldenrod'
  ) +
  facet_wrap(~cellID,scales = 'free')+
  theme_bw()+
  theme(axis.text.x = element_blank())
