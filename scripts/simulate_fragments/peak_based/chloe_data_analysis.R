library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
cell_scores = read.table("/data/rds/DMP/UCEC/GENEVOD/ccolson/MECCA/AllCells/HallmarkPathwayAnalysis/HIPEC_P05/C4_Normalised/magic_gsva_bootstrap_zscoreParam_consensus_pathway_scores.csv",header = T,sep=",")
cell_epigenetic_df = read.table("/data/rds/DMP/UCEC/GENEVOD/ccolson/MECCA/AllCells/HallmarkPathwayAnalysis/HIPEC_P05/C4_Normalised/magic_gsva_bootstrap_zscoreParam_consensus_clusters.csv",header = T,sep=",")

genes_hallmark_list = readRDS("/data/rds/DMP/UCEC/GENEVOD/ccolson/MECCA/AllCells/HallmarkPathwayAnalysis/HallmarkPathways.rds")
genes_hallmark <- data.frame(
  pathway = rep(names(genes_hallmark_list), lengths(genes_hallmark_list)),
  gene = unlist(genes_hallmark_list, use.names = FALSE)
)

peaks = readRDS("/data/rds/DMP/UCEC/GENEVOD/ccolson/MECCA/AllCells/ArchRProjects/HIPEC_P05_AllCells/PeakCalls/Type/Tumour-reproduciblePeaks.gr.rds")%>% as.data.frame()
peaks_df = peaks %>% 
  select(seqnames,start,end,nearestGene,peakType) %>% 
  dplyr::rename(gene=nearestGene) %>% 
  left_join(genes_hallmark,relationship = "many-to-many") %>% 
  filter(!is.na(pathway))

cell_scores_long = reshape2::melt(cell_scores)
colnames(cell_scores_long) = c("cell_id","pathway","score")


hallmark_map <- list(
  immune = c(
    "HALLMARK_ALLOGRAFT_REJECTION",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_IL2_STAT5_SIGNALING",
    "HALLMARK_IL6_JAK_STAT3_SIGNALING",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
  ),
  
  stromal = c(
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_ANGIOGENESIS",
    "HALLMARK_COMPLEMENT",
    "HALLMARK_COAGULATION",
    "HALLMARK_TGF_BETA_SIGNALING",
    "HALLMARK_HYPOXIA"
  ),
  
  oncogenic = c(
    "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
    "HALLMARK_MTORC1_SIGNALING",
    "HALLMARK_MYC_TARGETS_V1",
    "HALLMARK_MYC_TARGETS_V2",
    "HALLMARK_E2F_TARGETS",
    "HALLMARK_G2M_CHECKPOINT",
    "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
    "HALLMARK_KRAS_SIGNALING_UP",
    "HALLMARK_KRAS_SIGNALING_DN"
  ),
  
  stress = c(
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    "HALLMARK_DNA_REPAIR",
    "HALLMARK_P53_PATHWAY",
    "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
    "HALLMARK_REACTIVE_OXIGEN_SPECIES_PATHWAY"
  )
)

pathway_map_df <- enframe(hallmark_map, "class", "pathway") %>%
  unnest(pathway)

cell_scores_long =cell_scores_long %>% 
  left_join(pathway_map_df) %>% 
  mutate(class=case_when(is.na(class)~"other",
                         TRUE~class))

final_data =cell_scores_long %>% 
  left_join(y = cell_epigenetic_df,by="cell_id")
top_pathways =final_data %>% 
  group_by(pathway, consensus_cluster) %>% 
  summarise(
    mean_score = mean(score),
    sd_score = sd(score),
    .groups = "drop"
  ) %>% 
  group_by(pathway) %>% 
  filter(any(abs(mean_score) > 0.8))  %>% 
  pull(pathway) %>% unique()

final_data_filtered = final_data %>% 
  filter(pathway%in%top_pathways)


final_data_filtered_rescaled =final_data_filtered %>% 
  # group_by(pathway,consensus_cluster) %>%
  # summarise(mean_score=mean(score)) %>% 
  mutate(
    score_01 = scales::rescale(score, to = c(0, 1))
  ) #%>%
  # ungroup()

final_data_filtered_rescaled_per_group=final_data_filtered %>% 
  group_by(pathway,consensus_cluster) %>%
  summarise(mean_score=mean(score)) %>%
  mutate(
    mean_score_01 = scales::rescale(mean_score, to = c(0, 1))
  ) #%>%
# ungroup()


final_data_pathway =final_data_filtered_rescaled %>% 
  group_by(pathway,consensus_cluster) %>% 
  summarise(mean_score=mean(score),
            mean_score_01=mean(score_01))


mat_zcore <- final_data_pathway %>%
  select(consensus_cluster, pathway, mean_score) %>%
  pivot_wider(names_from = consensus_cluster, values_from = mean_score) %>% 
  column_to_rownames(var = 'pathway') %>% 
  as.matrix()


Heatmap(
  mat_zcore,
  name = "Score",
  show_row_names = TRUE,
  show_column_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = T,
)

mat_01 <- final_data_pathway %>%
  select(consensus_cluster, pathway, mean_score_01) %>%
  pivot_wider(names_from = consensus_cluster, values_from = mean_score_01) %>% 
  column_to_rownames(var = 'pathway') %>% 
  as.matrix()


Heatmap(
  mat_01,
  name = "Score",
  show_row_names = TRUE,
  show_column_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = T,
  col = colorRampPalette(c("white", "darkgreen"))(10)
)


mat_01_group <- final_data_filtered_rescaled_per_group %>%
  select(consensus_cluster, pathway, mean_score_01) %>%
  pivot_wider(names_from = consensus_cluster, values_from = mean_score_01) %>% 
  column_to_rownames(var = 'pathway') %>% 
  as.matrix()


Heatmap(
  mat_01_group,
  name = "Score",
  show_row_names = TRUE,
  show_column_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = T,
  col = colorRampPalette(c("white", "darkgreen"))(10)
)


dir.create(path = 'input_data_P05')
saveRDS(object = mat_01_group,file = 'input_data_P05/a_scores_scaled_per_group.rds')
saveRDS(object = mat_01,file = 'input_data_P05/a_scores_all.rds')


filtered_peaks = peaks_df %>% 
  filter(pathway%in%top_pathways)
changing_peaks = peaks_df %>% 
  filter(!pathway%in%top_pathways)
saveRDS(object = filtered_peaks,file = 'input_data_P05/peak_per_pathways.rds')



saveRDS(object = changing_peaks,file = 'input_data_P05/changing_peaks.rds')


