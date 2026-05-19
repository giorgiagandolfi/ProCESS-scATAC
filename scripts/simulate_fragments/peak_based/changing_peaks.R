rm(list=ls())
library(ProCESS)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC/scripts/plot_muller.R")
setwd("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-scATAC/scripts/simulate_fragments/peak_based/utils.R")
sample_forest <- load_sample_forest("sample_forest_atac_small.sff")
phylo_forest <- load_phylogenetic_forest("phylo_forest_atac_small.sff")


mutant_cols <- c(
  "Clone 1" = "goldenrod",
  "Clone 2" = "magenta4",
  "Clone 3" ="forestgreen"
)

sample_cols <- c(
  S1 = "darkorange",
  S2 = "royalblue3"
)

fixed_peaks <- paste0("peak_",seq_along(1:100))
status_fixed_peaks <- rbinom(n = 100, size = 1, prob = 0.4)
status_fixed_peaks <- rep(1,100)
df_fixed <- data.frame(peak_id=fixed_peaks,
                       status=status_fixed_peaks) 

####### define the markov chain of open-closed chromatin
total_simulated_peaks <- 10
fixed_peaks <- rep(x = "fixed",total_simulated_peaks*0.8)
changing_peaks <- rep(x = "fluctuating",total_simulated_peaks*0.2)

# df_peak_types <- data.frame(peak_id=paste0("peak_",seq_along(1:total_simulated_peaks)),
#                             type=c(fixed_peaks,changing_peaks)) 
crc_peaks<- readRDS("crc_peaks.rds")
df_peak_types <- crc_peaks %>% head(10) %>% 
  mutate(type=c(fixed_peaks,changing_peaks))
df_peak_types <- df_peak_types %>% 
  mutate(status=case_when(type=="fixed"~"open",
                          TRUE ~NA))


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
  print(cell_id_node)
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
      mutate(cell_id=cell_id_node)
    list_fragment_counts_per_peak <- list()
    for (p in 1:nrow(df_peak_types_cell)){
      peak_df = df_peak_types_cell[p,]
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

start <- Sys.time()
tour <- get_label_tour(sample_forest, labelling_functor1, only_leaves=TRUE)
end <- Sys.time()
end - start

df <- list()
i=1
start <- Sys.time()
while (!tour$done) {
  # print(tour$value)
  df[[i]] <- tour$value
  i=i+1
  tour$step()
}
end <- Sys.time()
end - start

df_final <- do.call("rbind",df)
cell_info <- sample_forest$get_nodes()
df_final <- df_final %>% inner_join(cell_info) %>% 
  select(!label.cell_id)
colnames(df_final) <- gsub("label.", "", colnames(df_final))
df_final<- df_final %>% 
  mutate(frag_final=paste0(fragment,"_",peak)) %>% 
  mutate(frag_size=to-from)




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
    data = df_final,
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
    data = selected_frags %>% filter(fragment_len <= 800),
    aes(x = fragment_len),
    binwidth = 1,
    # fill = "tomato",
    color="tomato",
    alpha = 0.3
  )+theme_minimal()


df_final %>% 
  group_by(peak,cell_id) %>% 
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

####################################

mutant_cols <- c(
  "A" = "goldenrod",
  "B" = "magenta4"
  # "Clone 3" ="forestgreen"
)

sample_cols <- c(
  S1 = "darkorange",
  S2 = "royalblue3"
)
df_final <- df_final %>% 
  # filter(!str_starts(label.peak_id, "peak_")) %>% 
  dplyr::select(cell_id,peak_id,status,mutant,sample) 
row_annot_df <- df_final %>%
  distinct(cell_id, sample, mutant) %>%
  arrange(sample)

df_final <- df_final %>% 
  mutate(status=case_when(status=="open"~1,
                                TRUE~0)) %>% 
  mutate(cell_id = factor(cell_id, levels = row_annot_df$cell_id)) %>%
  arrange(cell_id)
# convert to matrix: rows = cell_id, cols = peak_id
mat <- df_final %>%
  select(cell_id, peak_id, status) %>%
  unique() %>% 
  pivot_wider(
    names_from = peak_id,
    values_from = status
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
  Clone = row_annot_df$mutant,
  # sample =row_annot_df$sample,
  col = list(Clone=mutant_cols)
)
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
  left_annotation = row_ha,show_row_dend = F,
  show_column_names = F, show_row_names = F, split = row_annot_df$sample
)
# png(filename = "mets_plot/plot_heatmap_chromatin.png",
#     width = 6, height = 6, units = "in", res = 300)
draw(plot_heatmap_chromtin)
# dev.off()

plot_forest(sample_forest,color_map = mutant_cols)

