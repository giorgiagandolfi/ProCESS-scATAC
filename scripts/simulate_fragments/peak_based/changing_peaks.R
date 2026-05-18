rm(list=ls())
library(ProCESS)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC/scripts/plot_muller.R")
setwd("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC")
sample_forest <- load_sample_forest("sample_forest_atac3.sff")
phylo_forest <- load_phylogenetic_forest("phylo_forest_atac3.sff")


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
states <- c("open", "closed")

# Transition matrix:
# rows = current state, cols = next state
P <- matrix(c(
  0.8, 0.2,  # open -> open/closed
  0.3, 0.7   # closed -> open/closed
), nrow = 2, byrow = TRUE)

rownames(P) <- states
colnames(P) <- states
df_fixed <- data.frame(peak_id=paste0("peak_",seq_along(1:1000)),
                       cell_id=NA,
                       status=NA) 

simulate_fluctating_peaks <- function(n = 1000, P, states) {
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

labelling_functor1 <- function(label, node) {
  # the nodes are labelled by the identifiers of the associated cells
  cell_id_node = node$cell_id
  chromatin_status <- simulate_fluctating_peaks(n = 1000,P,states = c("open","closed"))
  df_changing <- data.frame(peak_id=paste0("peak_",seq_along(1:1000)),
                            cell_id=cell_id_node,
                            status=chromatin_status)
  # df_final <- rbind(df_fixed,df_changing)
  return(df_changing)  
}

tour <- get_label_tour(sample_forest, labelling_functor1, only_leaves=TRUE)
df <- list()
i=1
while (!tour$done) {
  # print(tour$value)
  df[[i]] <- tour$value
  i=i+1
  tour$step()
}
df_final <- do.call("rbind",df)
cell_info <- sample_forest$get_nodes()
df_final <- df_final %>% inner_join(cell_info)
df<-df_final





########## plot heatmap ###########



mutant_cols <- c(
  "Clone 1" = "goldenrod",
  "Clone 2" = "magenta4",
  "Clone 3" ="forestgreen"
)

sample_cols <- c(
  S1 = "darkorange",
  S2 = "royalblue3"
)
final_df <- df %>% 
  # filter(!str_starts(label.peak_id, "peak_")) %>% 
  select(cell_id,label.peak_id,label.status,mutant,sample) 
row_annot_df <- final_df %>%
  distinct(cell_id, sample, mutant) %>%
  arrange(sample)

final_df <- final_df %>% 
  mutate(label.status=case_when(label.status=="open"~1,
                                TRUE~0)) %>% 
  mutate(cell_id = factor(cell_id, levels = row_annot_df$cell_id)) %>%
  arrange(cell_id)
# convert to matrix: rows = cell_id, cols = peak_id
mat <- final_df %>%
  select(cell_id, label.peak_id, label.status) %>%
  pivot_wider(
    names_from = label.peak_id,
    values_from = label.status
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
  cluster_rows = T,
  cluster_columns = F,
  left_annotation = row_ha,show_row_dend = F,
  show_column_names = F, show_row_names = F, split = row_annot_df$sample
)
# png(filename = "mets_plot/plot_heatmap_chromatin.png",
#     width = 6, height = 6, units = "in", res = 300)
draw(plot_heatmap_chromtin)
# dev.off()

