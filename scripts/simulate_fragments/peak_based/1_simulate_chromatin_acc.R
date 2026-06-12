# rm(list=ls())
library(ProCESS)
library(dplyr)
# library(ggplot2)
library(tidyverse)
library(readxl)
# library(ComplexHeatmap)
# library(circlize)
# library(data.table)
#source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC/scripts/plot_muller.R")
#setwd("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-scATAC/scripts/simulate_fragments/peak_based/")
#source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-scATAC/scripts/simulate_fragments/peak_based/utils.R")
source("utils.R")

sample_forest <- load_sample_forest("sample_forest_atac_epigenome.sff")
phylo_forest <- load_phylogenetic_forest("phylo_forest_atac_epigenetic.sff")


mutant_cols <- c(
  "A+" = "goldenrod",
  "A-" = "magenta4",
  "B+" = "royalblue3",
  "B-" ="forestgreen"
)

sample_cols <- c(
  S1 = "lightskyblue1",
  S2 = "royalblue3"
)



####### define the markov chain of open-closed chromatin
# crc_peaks<- readRDS("crc_peaks.rds")
crc_peaks <- read_excel('/data/rds/DMP/UCEC/GENEVOD/ggandolfi/HTAN_data/41586_2023_6682_MOESM4_ESM.xlsx',sheet = 1) %>% 
  dplyr::filter(Cancer=="CRC") %>% 
  dplyr::filter(!grepl("chrX", peak)) %>% 
  dplyr::filter(!grepl("chrY", peak))
crc_peaks_top_fch <- read_excel('/data/rds/DMP/UCEC/GENEVOD/ggandolfi/HTAN_data/41586_2023_6682_MOESM4_ESM.xlsx',sheet = 3) %>% 
  filter(Cancer=="CRC") %>% 
  filter(!grepl("chrX", peak)) %>% 
  filter(!grepl("chrY", peak))
crc_peaks_not_fch = anti_join(crc_peaks,crc_peaks_top_fch)

total_simulated_peaks <- nrow(crc_peaks)
fixed_peaks <- rep(x = "fixed",nrow(crc_peaks_top_fch))
changing_peaks <- rep(x = "fluctuating",nrow(crc_peaks_not_fch)*0.4)
clone_peaks<- rep(x = "clonal",nrow(crc_peaks_not_fch)*0.6)

# fixed_peaks <- rep(x = "fixed",total_simulated_peaks*0.3)
# changing_peaks <- rep(x = "fluctuating",total_simulated_peaks*0.2)
# clone_peaks<- rep(x = "clonal",total_simulated_peaks*0.5)
# df_peak_types <- data.frame(peak_id=paste0("peak_",seq_along(1:total_simulated_peaks)),
#                             type=c(fixed_peaks,changing_peaks))

df_peak_types <- crc_peaks %>% 
  mutate(type=c(fixed_peaks,changing_peaks,clone_peaks))

df_peak_types <- df_peak_types %>% 
  group_by(type) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  mutate(mutant = NA_character_,
         epistate = NA_character_)
  
clonal_idx <- df_peak_types$type == "clonal"
nC <- sum(clonal_idx)

A <- round(nC * 0.6)
B <- nC - A

df_peak_types$mutant[clonal_idx] <- sample(c(rep("A", A), rep("B", B)))
df_peak_types$epistate[clonal_idx] <- sample(c("+", "-"), nC, replace = TRUE)
df_peak_types <- df_peak_types %>%
  mutate(status=case_when(type=="fixed"~"open",
                          TRUE ~"closed")) %>%
  separate(peak,into = c("chr","from","to"),sep = "-",remove = F) %>% 
  mutate(from=as.numeric(from)) %>% 
  mutate(to=as.numeric(to))

# df_peak_types <- df_peak_types %>% 
#   mutate(status=case_when(type=="fixed"~"open",
#                           TRUE ~"closed")) %>% 
#   mutate(mutant=c(rep(NA,100*0.5),rep("A",100*0.2),rep("B",100*0.3))) %>% 
#   mutate(epistate=c(rep(NA,100*0.5),rep("+",100*0.1),rep("-",100*0.1),rep("+",100*0.2),rep("-",100*0.1)))


df_peak_types %>% 
  filter(type=="clonal") %>% 
  mutate(clone=paste0(mutant,epistate)) %>% 
  ggplot(aes(x=mutant,fill=clone)) +
  geom_bar()+
  scale_fill_manual(values=mutant_cols)+
  theme_minimal()+
  ggtitle(label = "Peak distribution per class")

# Transition matrix:
# rows = current state, cols = next state
states <- c("open", "closed")
P <- matrix(c(
  0.8, 0.2,  # open -> open/closed
  0.3, 0.7   # closed -> open/closed
), nrow = 2, byrow = TRUE)

rownames(P) <- states
colnames(P) <- states
# df_final_counts <- data.frame(peak_id=paste0("peak_",seq_along(1:total_simulated_peaks)),
#                        cell_id=NA,
#                        status=NA) 

simulate_fluctating_peaks <- function(n = 1, P, states) {
  x <- character(n)
  x[1] <- sample(states, 1)
  
  for (i in 2:n) {
    current <- x[i - 1]
    x[i] <- sample(
      states,
      size = 1,
      prob = P[current, ]
    )
  }
  x
}


cna_data_sim <- phylo_forest$get_cell_allelic_fragmentation()
cna_data_sim = cna_data_sim %>% 
  mutate(chr=paste0('chr',chr))
#### convert to absolute coordinates
# genome_coordinates <- CNAqc::gene_coordinates_GRCh38
# vfrom = genome_coordinates$from
# names(vfrom) = genome_coordinates$chr
# cna_data_sim = cna_data_sim %>% 
#   mutate(chr=paste0('chr',chr)) %>% 
#   mutate(from = begin + vfrom[chr],
#          to = end + vfrom[chr])
labelling_functor1 <- function(label, node) {
  
  
  # the nodes are labelled by the identifiers of the associated cells
  cell_id_node = node$cell_id
  cell_mutant = phylo_forest$get_nodes() %>% 
    filter(cell_id==cell_id_node) %>% 
    pull(mutant)
  cell_phenotype = phylo_forest$get_nodes() %>% 
    filter(cell_id==cell_id_node) %>% 
    pull(epistate)
  
  cell_cna = cna_data_sim %>% filter(cell_id==cell_id_node) %>% 
    select(chr,major,minor,from,to)
  
  if (nrow(cell_cna)!=0){
    ### update the status dataframe for each peak
    tot_fluctutating_peaks <- df_peak_types %>% 
      filter(type=="fluctuating") %>% nrow()
    chromatin_status_flc_peaks <- simulate_fluctating_peaks(n = tot_fluctutating_peaks,P,states = c("open","closed"))
    df_peak_types_cell <- df_peak_types %>%
      mutate(
        status = replace(
          status,
          type == "fluctuating",
          chromatin_status_flc_peaks
        )
      ) %>% 
      mutate(status=case_when(mutant==cell_mutant & epistate==cell_phenotype ~ "open",TRUE~status)) %>% 
      mutate(cell_id=cell_id_node)
    
    
    list_fragment_counts_per_peak <- list()
    df_peak_types_cell_open <- df_peak_types_cell %>% 
      filter(status=="open")
    for (p in 1:nrow(df_peak_types_cell_open)){
      peak_df = df_peak_types_cell_open[p,]
      peak_status = peak_df %>% pull(status)
      # print(p)
      peak_df = peak_df %>%
        left_join(cell_cna, by = "chr") %>%
        filter(
          from.x <= to.y,
          to.x   >= from.y
        ) %>%
        transmute(peak,peak_id,chr,peak_from = from.x,
                  peak_to= to.x,major, minor,status, cell_id, type)
      total_cna_peak=peak_df$major+peak_df$minor
      print(total_cna_peak)
      
      fragment_counts <-sample_fragments_for_peak(peak_id =  peak_df %>% pull(peak_id) %>% as.character(),
                                                  peak_from =  peak_df %>% pull(peak_from) %>% as.numeric(),
                                                  peak_to = peak_df %>% pull(peak_to)%>% as.numeric(),
                                                  fragment_len_dist = selected_frags_dist_len,
                                                  tot_cn = total_cna_peak) %>%
        mutate(cell_id=cell_id_node) %>%
        dplyr::rename(peak_id=peak)
      
      list_fragment_counts_per_peak[[p]] <- full_join(fragment_counts,peak_df,by=c("cell_id","peak_id"))
    } 
  } else {
    tot_fluctutating_peaks <- df_peak_types %>% 
      filter(type=="fluctuating") %>% nrow()
    chromatin_status_flc_peaks <- simulate_fluctating_peaks(n = tot_fluctutating_peaks,P,states = c("open","closed"))
    df_peak_types_cell <- df_peak_types %>%
      mutate(
        status = replace(
          status,
          type == "fluctuating",
          chromatin_status_flc_peaks
        )
      ) %>% 
      mutate(cell_id=cell_id_node)
    list_fragment_counts_per_peak <- list()
    for (p in 1:nrow(df_peak_types_cell)){
      list_fragment_counts_per_peak[[p]]<-NA
    }
  }
  df_final <- do.call("rbind",list_fragment_counts_per_peak)
  return(df_final)  
}

labelling_functor2 <- function(label, node) {
  
  
  # the nodes are labelled by the identifiers of the associated cells
  # cell_id_node =115001
  cell_id_node = node$cell_id
  cell_mutant = phylo_forest$get_nodes() %>% 
    filter(cell_id==cell_id_node) %>% 
    pull(mutant)
  cell_phenotype = phylo_forest$get_nodes() %>% 
    filter(cell_id==cell_id_node) %>% 
    pull(epistate)
  
  cell_cna = cna_data_sim %>% filter(cell_id==cell_id_node) %>% 
    dplyr::rename(from=begin) %>% 
    dplyr::rename(to=end) %>% 
    select(chr,major,minor,from,to)
  
  if (nrow(cell_cna)!=0){
    ### update the status dataframe for each peak
    tot_fluctutating_peaks <- df_peak_types %>% 
      filter(type=="fluctuating") %>% nrow()
    chromatin_status_flc_peaks <- simulate_fluctating_peaks(n = tot_fluctutating_peaks,P,states = c("open","closed"))
    df_peak_types_cell <- df_peak_types %>%
      mutate(
        status = replace(
          status,
          type == "fluctuating",
          chromatin_status_flc_peaks
        )
      ) %>% 
      mutate(status=case_when(mutant==cell_mutant & epistate==cell_phenotype ~ "open",TRUE~status)) %>% 
      mutate(cell_id=cell_id_node)
    # df_final <- do.call("rbind",list_fragment_counts_per_peak)
    return(df_peak_types_cell)  
  }
}

labelling_functor3 <- function(label, node) {
  
  
  # the nodes are labelled by the identifiers of the associated cells
  cell_id_node = node$cell_id
  cell_mutant = phylo_forest$get_nodes() %>% 
    filter(cell_id==cell_id_node) %>% 
    pull(mutant)
  cell_phenotype = phylo_forest$get_nodes() %>% 
    filter(cell_id==cell_id_node) %>% 
    pull(epistate)
  
  cell_cna = cna_data_sim %>% filter(cell_id==cell_id_node) %>% 
    dplyr::rename(from=begin) %>% 
    dplyr::rename(to=end) %>% 
    select(chr,major,minor,from,to)
  
  if (nrow(cell_cna)!=0){
    ### update the status dataframe for each peak
    df_peak_types_cell <- df_peak_types %>%
      mutate(cell_id=cell_id_node)
    peaks_dt <- as.data.table(df_peak_types_cell)
    cna_dt   <- as.data.table(cell_cna)
    setkey(peaks_dt, chr, from, to)
    setkey(cna_dt, chr, from, to)
    res =foverlaps(
      peaks_dt,
      cna_dt,
      by.x = c("chr", "from", "to"),
      by.y = c("chr", "from", "to"),
      type = "any",   # overlap allowed
      nomatch = NA
    )
    df_peak_with_cna <- as.data.frame(res) %>% 
      mutate(tot_cn=major+minor) %>% 
      dplyr::select(cell_id,peak,tot_cn)
  } else {
    df_peak_with_cna <- df_peak_types %>%
      mutate(cell_id=cell_id_node) %>% 
      mutate(tot_cn=NA)
  }
  return(df_peak_with_cna)
}

start <- Sys.time()
tour_peaks <- get_label_tour(sample_forest, labelling_functor2, only_leaves=TRUE)
tour_cna <- get_label_tour(sample_forest, labelling_functor3, only_leaves=TRUE)
end <- Sys.time()
end - start

list_peaks_final <- list()
i=1
start <- Sys.time()
while (!tour_peaks$done) {
  # print(tour$value)
  list_peaks_final[[i]] <- tour_peaks$value
  i=i+1
  print(i)
  tour_peaks$step()
}
end <- Sys.time()
end - start
df_peak_final <- do.call("rbind",list_peaks_final)

list_peak_cna_final <- list()
i=1
start <- Sys.time()
while (!tour_cna$done) {
  # print(tour$value)
  list_peak_cna_final[[i]] <- tour_cna$value
  i=i+1
  tour_cna$step()
}
end <- Sys.time()
end - start
df_peak_cna_final <- do.call("rbind",list_peak_cna_final)


########

df_peak_final <- df_peak_final %>% dplyr::select(cell_id,label.peak,label.status) %>% 
  mutate(status=case_when(label.status=="open"~1,
                          TRUE~0)) %>% 
  dplyr::rename(peak_id=label.peak) %>% 
  dplyr::select(cell_id,peak_id,status)
df_peak_cna_final <- df_peak_cna_final %>% dplyr::select(cell_id,label.peak,label.tot_cn) %>% 
  dplyr::rename(peak_id=label.peak) %>% 
  dplyr::rename(tot_cn=label.tot_cn) %>% 
  dplyr::select(cell_id,peak_id,tot_cn)
cell_info <- sample_forest$get_nodes()
saveRDS(object = df_peak_cna_final,file = "df_peak_cna_final_big.rds")
saveRDS(object = df_peak_final,file = "df_peak_final_big.rds")

saveRDS(object = cell_info,file = "cell_info_big.rds")
