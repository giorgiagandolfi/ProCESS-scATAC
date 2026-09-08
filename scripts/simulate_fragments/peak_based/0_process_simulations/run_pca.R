library(dplyr)
library(uwot)
library(ggplot2)
library(patchwork)
library(ProCESS)
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project/ProCESS-scATAC/scripts/simulate_fragments/peak_based/utils.R")
outdir = "/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project/analysis/results_case1"
cell_ids <- list.files(path = paste0(outdir,"/","process_simulation"))
peaks_sim_df <- lapply(cell_ids, function(cell){
  readRDS(file.path(outdir,"process_simulation",cell,paste0("cell_",cell,"_peak_accessibility.rds")))
  # print(cell)
}) %>% bind_rows() 
peaks_sim_df_open <- peaks_sim_df %>% 
  select(cell_id,peak,status)
peak_matrix <- xtabs(status ~ cell_id + peak, data = peaks_sim_df_open)
peak_matrix <- as.matrix(peak_matrix)
pca <- prcomp(peak_matrix, scale. = TRUE)


pca_df <- data.frame(
  cell_id = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
)
peaks_sim_df <- peaks_sim_df %>% 
  mutate(cell_id=as.character(cell_id))
pca_df <- pca_df %>% 
  inner_join(peaks_sim_df,by="cell_id")


color_map_epi = c("E1"="forestgreen",
                  "E2"="goldenrod",
                  "E3"="orchid2"
)
color_map_gen = c("G1"="coral2",
                  "G2"="turquoise4",
                  "G3"="darkorange")
p1=pca_df %>% 
  select(PC1,PC2,cell_id,mutant,epistate) %>% 
  distinct() %>% 
  ggplot(aes(x = PC1, y = PC2,color=mutant)) +
  geom_point(size = .5) +
  theme_minimal()+
  # theme_classic()+
  scale_color_manual(values = color_map_gen)+
  labs(
    title = "Cell coloured by genetic",
    x = "PC1",
    y = "PC2"
  )
p2=pca_df %>% 
  select(PC1,PC2,cell_id,mutant,epistate) %>% 
  distinct() %>% 
  ggplot(aes(x = PC1, y = PC2,color=epistate)) +
  geom_point(size = 0.5) +
  theme_minimal()+
  # theme_classic()+
  scale_color_manual(values = color_map_epi)+
  labs(
    title = "Cell coloured by epistate",
    x = "PC1",
    y = "PC2"
  )
p1+p2
simu_path <- readRDS("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project/reference_data/a_scores_pathway_case1.rds")
path_df=convert_activity_list(simu_path)

library(tidyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)

wide <- path_df %>%
  pivot_wider(names_from = epistate, values_from = activity) %>%
  column_to_rownames("pathway")

mat <- as.matrix(wide)
PEAK_FILE=readRDS("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project/reference_data/peak_pathway_list_unique_peaks.rds")

col_fun <- colorRamp2(c(0, 1), c("white", "#187fc4e1"))
Heatmap(t(mat), col = col_fun, name = "activity",show_column_names  = F,show_column_dend = F)


sample_forest <- load_sample_forest("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project/ProCESS-scATAC/scripts/simulate_fragments/peak_based/0_process_simulations/sample_forest_atac_case_1.sff")
color_map_epi = c("G1[E1]"="forestgreen",
                  "G2[E1]"="forestgreen",
                  "G3[E1]"="forestgreen",
                  "G1[E2]"="goldenrod",
                  "G2[E2]"="goldenrod",
                  "G3[E2]"="goldenrod",
                  "G1[E3]"="orchid2",
                  "G2[E3]"="orchid2",
                  "G3[E3]"="orchid2"
)
color_map_gen = c("G1[E1]"="coral2",
                  "G2[E1]"="turquoise4",
                  "G3[E1]"="darkorange",
                  "G1[E2]"="coral2",
                  "G2[E2]"="turquoise4",
                  "G3[E2]"="darkorange",
                  "G1[E3]"="coral2",
                  "G2[E3]"="turquoise4",
                  "G3[E3]"="darkorange"
)
p1_forest =plot_forest(sample_forest,color_map = color_map_gen)+ggtitle(label = "Cells coloured by genetic")
p2_forest=plot_forest(sample_forest,color_map = color_map_epi)+ggtitle(label = "Cells coloured by epistate")
p1_forest+p2_forest
