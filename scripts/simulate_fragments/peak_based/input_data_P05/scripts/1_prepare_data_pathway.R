library(dplyr)
library(tidyverse)
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

genes_involvement = genes_hallmark %>% dplyr::group_by(gene) %>% dplyr::summarise(n_pathways=dplyr::n())

pathways <- names(genes_hallmark_list)
overlap_matrix <- matrix(
  0,
  nrow = length(pathways),
  ncol = length(pathways),
  dimnames = list(pathways, pathways)
)

for (x in pathways){
  for (y in pathways){
    common_genes <-length(
      intersect(genes_hallmark_list[[x]], genes_hallmark_list[[y]])
    )
    overlap_matrix[x, y] <- common_genes/length(genes_hallmark_list[[x]])
  }
}
Heatmap(matrix = overlap_matrix,col = colorRampPalette(c("white", "darkgreen"))(10),cluster_columns = F,cluster_rows = F)




epigenetic_scores = cell_scores %>% 
  dplyr::rename(cell_id=X) %>% 
  pivot_longer(
    cols = -cell_id,
    names_to = "pathway",
    values_to = "score"
  ) %>% 
  dplyr::inner_join(cell_epigenetic_df) %>% 
  dplyr::group_by(consensus_cluster,pathway) %>% 
  dplyr::summarise(mean_score=mean(score),
                   sd_score=sd(score)) %>% 
  dplyr::group_by(pathway) %>% 
  dplyr::mutate(
    mean_score_01_norm = (mean_score - min(mean_score)) /
      (max(mean_score) - min(mean_score))
  ) %>% 
  dplyr::mutate(
    mean_score_01_nor_offset = 0.05 + 
      ((mean_score - min(mean_score)) /
         (max(mean_score) - min(mean_score))) * 0.95
  ) %>% 
  dplyr::mutate(mean_score_logist=plogis(mean_score))



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

epigenetic_scores =epigenetic_scores %>% 
  left_join(pathway_map_df) %>% 
  mutate(class=case_when(is.na(class)~"other",
                         TRUE~class))

# top_pathways =epigenetic_scores %>% 
#   filter(sd_score>1)  %>% 
#   pull(pathway) %>% unique()


epigenetic_cluster_cols = c(
  "1" = "olivedrab",
  "2" = "forestgreen",
  "3" = "darkgreen",
  "4" ="palegreen1"
)

mat_zcore <- epigenetic_scores %>%
  select(consensus_cluster, pathway, mean_score) %>%
  pivot_wider(names_from = consensus_cluster, values_from = mean_score) %>% 
  column_to_rownames(var = 'pathway') %>% 
  as.matrix()

col_annot = HeatmapAnnotation(df = epigenetic_scores %>% ungroup() %>% select(consensus_cluster) %>% distinct(),col = list(
  consensus_cluster = epigenetic_cluster_cols
))
annotation <- epigenetic_scores %>% 
  select(pathway,class) %>% 
  distinct()

annotation <- annotation[match(rownames(mat_zcore), annotation$pathway), ]

row_ann = rowAnnotation(class = annotation$class)

mat_zcore_ht=Heatmap(
  mat_zcore,
  name = "Z-Score",
  show_row_names = TRUE,
  show_column_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = T,
  top_annotation = col_annot,left_annotation = row_ann
)

mat_01core <- epigenetic_scores %>%
  select(consensus_cluster, pathway, mean_score_01_norm) %>%
  pivot_wider(names_from = consensus_cluster, values_from = mean_score_01_norm) %>% 
  column_to_rownames(var = 'pathway') %>% 
  as.matrix()

mat_01core_ht=Heatmap(
  mat_01core,
  name = "01 Score",
  show_row_names = TRUE,
  show_column_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = T,
  top_annotation = col_annot,col = colorRampPalette(c("white", "darkgreen"))(10)
)


mat_logit <- epigenetic_scores %>%
  select(consensus_cluster, pathway, mean_score_logist) %>%
  pivot_wider(names_from = consensus_cluster, values_from = mean_score_logist) %>% 
  column_to_rownames(var = 'pathway') %>% 
  as.matrix()

mat_logit_ht=Heatmap(
  mat_logit,
  name = "01 Score",
  show_row_names = TRUE,
  show_column_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = T,
  top_annotation = col_annot,col = colorRampPalette(c("white", "darkgreen"))(10)
)

pdf("merged_zscore_and_rescaled_heatmaps.pdf",width = 20,height = 15)
draw(mat_zcore_ht+mat_01core_ht+mat_logit_ht,heatmap_legend_side = "bottom",annotation_legend_side = "bottom",
     merge_legends = T)
dev.off()
activity_df=epigenetic_scores %>%
  select(consensus_cluster, pathway, mean_score_logist) %>% 
  mutate(epistate=case_when(consensus_cluster==1~"E1",
                            consensus_cluster==3~"E3",
                            consensus_cluster==2~"E2",
                            consensus_cluster==4~"E4",
                            TRUE~NA)) %>% 
  dplyr::rename(a_score=mean_score_logist) %>% 
  dplyr::select(pathway,a_score,epistate)


epistates <- activity_df$epistate %>% unique()

activity_list <- list()
for (epi in epistates){
    pathways_state =activity_df %>% 
      filter(epistate==epi) 
    pathways_vect <- pathways_state$a_score
    names(pathways_vect) <- pathways_state$pathway
    activity_list[[epi]]<-pathways_vect
}


saveRDS(object = activity_list,file = 'input_data_P05/data/final_version/a_scores_pathway.rds')




peaks = readRDS("/data/rds/DMP/UCEC/GENEVOD/ccolson/MECCA/AllCells/ArchRProjects/HIPEC_P05_AllCells/PeakCalls/Type/Tumour-reproduciblePeaks.gr.rds")%>% as.data.frame()
selected_peaks = peaks %>% 
  filter(distToTSS<=1000) %>% ### do not know if to use the score or the distance from the gene start
  # select(seqnames,start,end,nearestGene,peakType) %>% 
  dplyr::rename(gene=nearestGene)

unique_sim_peaks = peaks_df %>% pull(peak) %>% unique() %>% length()

# peaks_df %>% 
#   group_by(peak) %>% 
#   mutate(n_involved_pathways=n()) %>% 
#   mutate(n_involved_pathways = case_when(
#     n_involved_pathways >= 7 ~ ">=7",
#     TRUE ~ as.character(n_involved_pathways)
#   )) %>% 
#   group_by(n_involved_pathways) %>% 
#   summarise(pct=n()/unique_sim_peaks) %>% 
#   ggplot(aes(x = "", y = pct, fill = as.factor(n_involved_pathways))) +
#   geom_col(width = 1) +
#   coord_polar(theta = "y") +
#   theme_void() +
#   scale_fill_brewer(palette = 'Dark2')+
#   labs(fill = "Involved pathways")+
#   theme(legend.position = "bottom")


peaks_df_gene_filtered = peaks %>% 
  filter(distToTSS<=1000) %>% ### do not know if to use the score or the distance from the gene start
  # select(seqnames,start,end,nearestGene,peakType) %>% 
  dplyr::rename(gene=nearestGene) %>% 
  dplyr::left_join(genes_hallmark,relationship = "many-to-many") %>% 
  dplyr::filter(!is.na(pathway)) %>% 
  dplyr::mutate(peak=paste(seqnames,start,end,sep='-')) %>% 
  left_join(y = genes_involvement) %>% 
  left_join(activity_df,relationship = "many-to-many") %>% 
  dplyr::select(peak,gene,pathway,a_score,n_pathways,epistate)
  
  
peaks_df_gene_common = peaks_df_gene_filtered %>% 
  dplyr::filter(n_pathways>1) %>% 
  dplyr::group_by(pathway,gene) %>% 
  dplyr::summarise(mean_ascore=mean(a_score),
                   sd_ascore=sd(a_score)) %>% 
  group_by(gene) %>% 
  slice_max(order_by = sd_ascore,n = 1) %>% 
  inner_join(y = selected_peaks) %>% 
  dplyr::mutate(peak=paste(seqnames,start,end,sep='-')) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(peak,pathway)

peaks_df_gene_single = peaks_df_gene_filtered %>% 
  dplyr::filter(n_pathways==1) %>% 
  dplyr::select(peak,pathway,gene) %>% 
  distinct() %>% 
  left_join(y = selected_peaks,relationship = "many-to-many") %>% 
  dplyr::mutate(peak=paste(seqnames,start,end,sep='-')) %>% 
  dplyr::select(peak,pathway) %>% 
  distinct()

peaks_df = rbind(peaks_df_gene_single,peaks_df_gene_common)

pathways = peaks_df$pathway %>% unique()
peak_pathway_list = list()
for (pat in pathways){
  peak_pathway_list[[pat]] <- peaks_df %>% 
    filter(pathway==pat) %>% 
    separate(peak,into = c("chr","from","to"),sep = "-",remove = F,convert = T) %>% 
    mutate(peak_lenght=to-from) %>% 
    mutate(chr = str_remove(chr, "chr"))
}
saveRDS(object = peak_pathway_list,file = 'input_data_P05/data/final_version/peak_pathway_list_unique_peaks.rds')
