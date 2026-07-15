library(dplyr)
library(ComplexHeatmap)
library(circlize)
# library(ArchR)
library(tidyverse)
# 
# set.seed(1)
# addArchRThreads(threads = 16) 
# addArchRGenome("hg38")
cell_scores = read.table("/data/rds/DMP/UCEC/GENEVOD/ccolson/MECCA/AllCells/HallmarkPathwayAnalysis/HIPEC_P05/C4_Normalised/magic_gsva_bootstrap_zscoreParam_consensus_pathway_scores.csv",header = T,sep=",")
cell_epigenetic_df = read.table("/data/rds/DMP/UCEC/GENEVOD/ccolson/MECCA/AllCells/HallmarkPathwayAnalysis/HIPEC_P05/C4_Normalised/magic_gsva_bootstrap_zscoreParam_consensus_clusters.csv",header = T,sep=",")

genes_hallmark_list = readRDS("/data/rds/DMP/UCEC/GENEVOD/ccolson/MECCA/AllCells/HallmarkPathwayAnalysis/HallmarkPathways.rds")
genes_hallmark <- data.frame(
  pathway = rep(names(genes_hallmark_list), lengths(genes_hallmark_list)),
  gene = unlist(genes_hallmark_list, use.names = FALSE)
)
normalised_gene_score <- read.csv("/data/rds/DMP/UCEC/GENEVOD/ccolson/MECCA/AllCells/HallmarkPathwayAnalysis/HIPEC_P05/AllTumourCells_Gene_Normalised_Zscores_10-03-26.csv",row.names = 1)

mat_long <- normalised_gene_score %>% 
  rownames_to_column(var = "cell_id")  %>% 
  pivot_longer(
    cols = -cell_id,
    names_to = "gene",
    values_to = "score"
  ) %>% 
  mutate(score=as.numeric(score))

peaks = readRDS("/data/rds/DMP/UCEC/GENEVOD/ccolson/MECCA/AllCells/ArchRProjects/HIPEC_P05_AllCells/PeakCalls/Type/Tumour-reproduciblePeaks.gr.rds")%>% as.data.frame()
peaks_df = peaks %>% 
  filter(distToTSS<=1000) %>% ### do not know if to use the score or the distance from the gene start
  select(seqnames,start,end,nearestGene,peakType) %>% 
  dplyr::rename(gene=nearestGene) %>% 
  left_join(genes_hallmark,relationship = "many-to-many") %>% 
  filter(!is.na(pathway))

# peaks_df_with_gene_scores <- peaks_df %>% 
#   left_join(mat_long,relationship = "many-to-many") %>% 
#   left_join(cell_epigenetic_df)
# 
# peaks_df_with_gene_scores <- peaks_df_with_gene_scores %>% 
#   filter(!is.na(consensus_cluster))
# 
# mean_scores <- peaks_df_with_gene_scores %>% 
#   # left_join(cell_epigenetic_df,relationship = "many-to-many") %>% 
#   filter(!is.na(consensus_cluster)) %>% 
#   group_by(gene,consensus_cluster) %>% 
#   summarise(mean_score=mean(score)) %>% 
#   group_by(gene) %>%
#   mutate(
#     mean_score_norm = (mean_score - min(mean_score)) /
#       (max(mean_score) - min(mean_score))
#     # mean_score_z = (mean_score - mean(mean_score)) /
#     #   sd(mean_score)
#   )

mean_scores <- mat_long %>% 
  dplyr::left_join(cell_epigenetic_df,relationship = "many-to-many") %>%
  dplyr::filter(gene%in%peaks_df$gene) %>% 
  dplyr::filter(gene!="CBS") %>% 
  dplyr::filter(!is.na(consensus_cluster)) %>%
  dplyr::group_by(gene,consensus_cluster) %>%
  dplyr::summarise(mean_score=mean(score),
                   log1p_score=plogis(scale(log1p(score)))) %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(mean_score_norm =(mean_score - min(mean_score)) /
                  (max(mean_score) - min(mean_score))) %>% 
  dplyr::group_by(gene) %>%
  dplyr::mutate(
      mean_score_norm_per_epi = (mean_score - min(mean_score)) /
        (max(mean_score) - min(mean_score))
    )


mean_scores_logscaled <- mat_long %>% 
  dplyr::left_join(cell_epigenetic_df,relationship = "many-to-many") %>%
  dplyr::filter(gene%in%peaks_df$gene) %>% 
  dplyr::filter(gene!="CBS") %>% 
  dplyr::filter(!is.na(consensus_cluster)) %>%
  dplyr::group_by(gene,consensus_cluster) %>%
  dplyr::summarise(mean_score=mean(score)) %>% 
  ungroup() %>% 
  dplyr::mutate(log1p_score_sigmoid=plogis(scale(log1p(mean_score)))) %>% 
  dplyr::mutate(log1p_score=scales::rescale(log1p(mean_score), to = c(0, 1))) %>% 
  dplyr::mutate(
    mean_score_norm_per_epi = (mean_score - min(mean_score)) /
      (max(mean_score) - min(mean_score))
  )

scores_peaks_df = mean_scores_logscaled %>% 
  select(gene,consensus_cluster,log1p_score) %>% 
  left_join(peaks_df,relationship = "many-to-many")

saveRDS(object = scores_peaks_df,file = "../data/a_scores_per_gene_log1p_score.rds")
  
mat_mean_scores_norm_group <- mean_scores_logscaled %>%
  ungroup() %>%
  select(gene,mean_score_norm_per_epi,consensus_cluster) %>% 
  distinct() %>% 
  tidyr::pivot_wider(
    names_from = consensus_cluster,
    values_from = mean_score_norm_per_epi,values_fill = 0
  ) %>%
  as.data.frame()

mat_mean_scores_log1 <- mean_scores_logscaled %>%
  ungroup() %>%
  select(gene,log1p_score,consensus_cluster) %>% 
  distinct() %>% 
  tidyr::pivot_wider(
    names_from = consensus_cluster,
    values_from = log1p_score,values_fill = 0
  ) %>%
  as.data.frame()


mat_mean_scores_not_norm <- mean_scores_logscaled %>%
  ungroup() %>%
  select(gene,mean_score,consensus_cluster) %>% 
  distinct() %>% 
  tidyr::pivot_wider(
    names_from = consensus_cluster,
    values_from = mean_score,values_fill = 0
  ) %>%
  as.data.frame()

# set rownames
rownames(mat_mean_scores_norm_group) <- mat_mean_scores_norm_group$gene
mat_mean_scores_norm_group$gene <- NULL
mat_mean_scores_norm_group <- as.matrix(mat_mean_scores_norm_group)

# set rownames
rownames(mat_mean_scores_log1) <- mat_mean_scores_log1$gene
mat_mean_scores_log1$gene <- NULL
mat_mean_scores_log1 <- as.matrix(mat_mean_scores_log1)

# set rownames
rownames(mat_mean_scores_not_norm) <- mat_mean_scores_not_norm$gene
mat_mean_scores_not_norm$gene <- NULL
mat_mean_scores_not_norm <- as.matrix(mat_mean_scores_not_norm)


sampled_genes <- c(sample(x = rownames(mat_mean_scores_not_norm),size = 100),"MYC")
col_ann <- HeatmapAnnotation(epigenetic_class=colnames(mat_mean_scores_not_norm),
                         col=list("1"='forestgreen',"2"='purple','3'='goldenrod','4'='navyblue'))


ht_log1p=Heatmap(
  mat_mean_scores_log1[sampled_genes,],
  name = "log1p mean score",
  cluster_rows = F,
  cluster_columns = F,
  show_row_names = F,
  show_column_names = T,
  # top_annotation = col_ann,
  col = colorRampPalette(c("white", "darkgreen"))(10)
)


ht_01_norm=Heatmap(
  mat_mean_scores_norm_group[sampled_genes,],
  name = "01 normalised mean score",
  cluster_rows = F,
  cluster_columns = F,
  show_row_names = F,
  show_column_names = T,
  col = colorRampPalette(c("white", "darkorange"))(10)
)


ht_gene_score=Heatmap(
  mat_mean_scores_not_norm[sampled_genes,],
  name = "mean_score",
  cluster_rows = F,
  cluster_columns = F,
  show_row_names = F,
  show_column_names = T,
  col = colorRampPalette(c("white", "navyblue"))(10)
)

ht_list=ht_gene_score + ht_01_norm + ht_log1p
draw(ht_list)

peaks_df_with_gene_scores %>% 
  select(gene,score,consensus_cluster) %>% 
  distinct() %>% 
  filter(gene=="MYC") %>% 
  ggplot(aes(x=as.factor(consensus_cluster),y=score,fill=as.factor(consensus_cluster)))+
  geom_boxplot()

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
# old method
# top_pathways =final_data %>% 
#   group_by(pathway, consensus_cluster) %>% 
#   summarise(
#     mean_score = mean(score),
#     sd_score = sd(score),
#     .groups = "drop"
#   ) %>% 
#   group_by(pathway) %>% 
#   filter(any(abs(mean_score) > 0.8))  %>% 
#   pull(pathway) %>% unique()


top_pathways =final_data %>% 
  group_by(pathway) %>% 
  summarise(
    mean_score = mean(score),
    sd_score = sd(score),
    .groups = "drop"
  ) %>% 
  filter(sd_score>1)  %>% 
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


mat_zcore_ht=Heatmap(
  mat_zcore,
  name = "Score",
  show_row_names = TRUE,
  show_column_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = T,
)
draw(mat_zcore_ht,heatmap_legend_side = "bottom")
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

annotation <-final_data_filtered_rescaled_per_group %>% 
  # select(pathway, mean_score_01) %>% 
  left_join(pathway_map_df) %>% 
  mutate(class=case_when(is.na(class)~"other",
                         TRUE~class)) %>% 
  select(pathway,class) %>% 
  distinct() 
  
annotation <- annotation[match(rownames(mat_01_group), annotation$pathway), ]

row_ann = rowAnnotation(class = annotation$class)
mat_01_group_ht =Heatmap(
  mat_01_group,
  name = "Score",
  show_row_names = TRUE,
  show_column_names = FALSE,
  cluster_rows = T,
  cluster_columns = T,
  col = colorRampPalette(c("white", "darkgreen"))(10),
  left_annotation = row_ann
)
draw(mat_01_group_ht,heatmap_legend_side = "bottom")

dir.create(path = 'input_data_P05')
saveRDS(object = mat_01_group,file = 'input_data_P05/a_scores_scaled_per_group_filtered_tss1kb_peaks.rds')
saveRDS(object = mat_01,file = 'input_data_P05/a_scores_all_filtered_tss1kb_peaks.rds')


filtered_peaks = peaks_df %>% 
  filter(pathway%in%top_pathways)
changing_peaks = peaks_df %>% 
  filter(!pathway%in%top_pathways)
saveRDS(object = filtered_peaks,file = 'input_data_P05/peak_per_pathways_filtered_tss1kb_peaks.rds')



saveRDS(object = changing_peaks,file = 'input_data_P05/changing_peaks.rds')

##### gene distribution across pathways
genes_hallmark %>% filter(pathway%in%top_pathways) %>% 
  group_by(gene) %>% summarise(n_pathway=n()) %>% 
  group_by(n_pathway) %>% summarise(n=n()) %>% 
  ggplot(aes(x=as.factor(n_pathway),y=n))+geom_col()+theme_minimal()


##### most spread genes in hallmarks
genes_hallmark %>% filter(pathway%in%top_pathways) %>% 
  group_by(gene) %>% mutate(n_pathway=n()) %>% 
  filter(n_pathway>=6) %>% 
  left_join(pathway_map_df) %>%
  ggplot(aes(y=gene,fill=class))+geom_bar()+
  theme_minimal()+
  theme(legend.position = 'bottom')
  
