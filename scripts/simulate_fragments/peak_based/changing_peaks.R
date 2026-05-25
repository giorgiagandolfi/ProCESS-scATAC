# rm(list=ls())
library(ProCESS)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(data.table)
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC/scripts/plot_muller.R")
setwd("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-scATAC/scripts/simulate_fragments/peak_based/")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-scATAC/scripts/simulate_fragments/peak_based/utils.R")
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
total_simulated_peaks <- 100
fixed_peaks <- rep(x = "fixed",total_simulated_peaks*0.3)
changing_peaks <- rep(x = "fluctuating",total_simulated_peaks*0.2)
clone_peaks<- rep(x = "clonal",total_simulated_peaks*0.5)

# df_peak_types <- data.frame(peak_id=paste0("peak_",seq_along(1:total_simulated_peaks)),
#                             type=c(fixed_peaks,changing_peaks)) 
crc_peaks<- readRDS("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC/crc_peaks.rds")
df_peak_types <- crc_peaks %>% head(100) %>% 
  mutate(type=c(fixed_peaks,changing_peaks,clone_peaks))
df_peak_types <- df_peak_types %>% 
  mutate(status=case_when(type=="fixed"~"open",
                          TRUE ~"closed")) %>% 
  mutate(mutant=c(rep(NA,100*0.5),rep("A",100*0.2),rep("B",100*0.3))) %>% 
  mutate(epistate=c(rep(NA,100*0.5),rep("+",100*0.1),rep("-",100*0.1),rep("+",100*0.2),rep("-",100*0.1)))
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

#### convert to absolute coordinates
genome_coordinates <- CNAqc::gene_coordinates_GRCh38
vfrom = genome_coordinates$from
names(vfrom) = genome_coordinates$chr
cna_data_sim = cna_data_sim %>% 
  mutate(chr=paste0('chr',chr)) %>% 
  mutate(from = begin + vfrom[chr],
         to = end + vfrom[chr])
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
  # cell_id_node =82903
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
      dplyr::select(cell_id,peak_id,tot_cn)
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

df_peak_final <- df_peak_final %>% dplyr::select(cell_id,label.peak_id,label.status) %>% 
  mutate(status=case_when(label.status=="open"~1,
                          TRUE~0)) %>% 
  dplyr::rename(peak_id=label.peak_id) %>% 
  dplyr::select(cell_id,peak_id,status)
df_peak_cna_final <- df_peak_cna_final %>% dplyr::select(cell_id,label.peak_id,label.tot_cn) %>% 
  dplyr::rename(peak_id=label.peak_id) %>% 
  dplyr::rename(tot_cn=label.tot_cn) %>% 
  dplyr::select(cell_id,peak_id,tot_cn)



frag_res_all_cells <- lapply(df_peak_final$cell_id %>% unique(), function(c){
  peak_cell <- df_peak_final %>% filter(cell_id==c) %>% 
    filter(status==1) %>% 
    separate(peak_id,into = c("chr","from","to"),sep = ":",remove = F)
  peak_cna_cell <- df_peak_cna_final %>% 
    inner_join(peak_cell)
  fragm_df_cell <- sample_fragments_for_peak_vec(
    peak_id   = peak_cell$peak_id,
    peak_from = peak_cell$from %>% as.numeric(),
    peak_to   = peak_cell$to %>% as.numeric(),
    fragment_len_dist = selected_frags_dist_len,
    tot_cn    = peak_cna_cell$tot_cn%>% as.numeric(),cell_id = c
  ) %>% as.data.frame() 
  return(fragm_df_cell)
}) %>% bind_rows()





cell_info <- sample_forest$get_nodes()
frag_res_all_cells <- frag_res_all_cells %>% inner_join(cell_info)
colnames(df_final) <- gsub("label.", "", colnames(df_final))
frag_res_all_cells<- frag_res_all_cells %>% 
  mutate(frag_id=paste0(fragment_start,":",fragment_end))

frag_res_all_cells_all_status <- frag_res_all_cells
# frag_res_all_cells_all_status <- df_peak_final %>% 
#   left_join(frag_res_all_cells) 
########## plot heatmap ###########



starting_dist_frag_sizes <- replicate(
  100000,
  sample_fragment_size()
)

library(dplyr)
library(ggplot2)

ggplot() +
  
  # simulated distribution
  geom_histogram(
    data = frag_res_all_cells_all_status,
    aes(x = fragment_size, y = after_stat(density)),
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
    data = selected_frags %>% filter(fragment_len <= 800),
    aes(x = fragment_len),
    binwidth = 1,
    # fill = "tomato",
    color="tomato",
    alpha = 0.3
  )+theme_minimal()




####################################

mutant_cols <- c(
  "A" = "goldenrod",
  "B" = "magenta4"
  # "Clone 3" ="forestgreen"
)

epistate_cols <- c(
  "+" = "forestgreen",
  "-" = "darkblue"
  # "Clone 3" ="forestgreen"
)

mutant_cols <- c(
  "A+" = "goldenrod",
  "A-" = "magenta4",
  "B+" = "royalblue3",
  "B-" ="forestgreen"
)
peak_class <- c(
  clonal = "violet",
  fluctuating = "royalblue3",
  fixed = "darkorange"
)
frag_res_all_cells_all_status <- frag_res_all_cells_all_status %>% 
  mutate(status=1) %>% 
  # filter(!str_starts(label.peak_id, "peak_")) %>% 
  dplyr::select(cell_id,peak_id,mutant,sample,epistate,status) 
row_annot_df <- frag_res_all_cells_all_status %>%
  mutate(clone=paste0(mutant,epistate)) %>% 
  distinct(cell_id, sample,clone) %>%
  arrange(sample)
# col_annot_df <- df_final %>% 
#   distinct(peak_id,type)
frag_res_all_cells_all_status <- frag_res_all_cells_all_status %>% 
  mutate(cell_id = factor(cell_id, levels = row_annot_df$cell_id)) %>%
  arrange(cell_id)
# convert to matrix: rows = cell_id, cols = peak_id
mat <- frag_res_all_cells_all_status %>%
  select(cell_id, peak_id, status) %>%
  unique() %>% 
  pivot_wider(
    names_from = peak_id,
    values_from = status,values_fill = 0
  ) %>%
  column_to_rownames("cell_id") %>%
  as.matrix()

# # row annotation dataframe
# row_annot_mut_df <- final_df %>%
#   distinct(cell_id, mutant)# %>%
# #  arrange(cell_id)
# 
# row_annot_sample_df <- final_df %>%
#   distinct(cell_id, sample) #%>%
#   # arrange(cell_id)


# colors for mutant annotation


# row annotation


row_ha <- rowAnnotation(
  Clone = row_annot_df$clone,
  # Epistate = row_annot_df$epistate,
  # sample =row_annot_df$sample,
  col = list(Clone=mutant_cols)
)

# col_ha <- columnAnnotation(
#   Peak_Class = col_annot_df$type,
#   col = list(Peak_Class=peak_class)
# )

# heatmap colors
col_fun <- c(
  "0" = "white",
  "1" = "grey"
  
)

# draw heatmap


plot_heatmap_chromtin <-Heatmap(
  mat,
  name = "Chromatin status",
  col = col_fun,
  cluster_rows = F,
  cluster_columns = F,
  # top_annotation = col_ha,
  left_annotation = row_ha,show_row_dend = F,
  show_column_names = F, show_row_names = F#, split = row_annot_df$sample
)
# png(filename = "mets_plot/plot_heatmap_chromatin.png",
#     width = 6, height = 6, units = "in", res = 300)
draw(plot_heatmap_chromtin)
# dev.off()


plot_forest(sample_forest,color_map = mutant_cols)

