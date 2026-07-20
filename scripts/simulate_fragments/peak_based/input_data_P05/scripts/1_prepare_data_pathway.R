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
  )



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

top_pathways =epigenetic_scores %>% 
  filter(sd_score>1)  %>% 
  pull(pathway) %>% unique()


final_data_filtered = final_data %>% 
  filter(pathway%in%top_pathways)

epigenetic_cluster_cols = c(
  "1" = "#ABB6FE",
  "2" = "#465efdff",
  "3" = "#0117a7ff",
  "4" ="#9402eeff"
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

draw(mat_zcore_ht+mat_01core_ht,heatmap_legend_side = "bottom",annotation_legend_side = "bottom",
     merge_legends = T)

activity_df=epigenetic_scores %>%
  select(consensus_cluster, pathway, mean_score_01_norm) %>% 
  mutate(epistate=case_when(consensus_cluster==1~"E1",
                            consensus_cluster==3~"E3",
                            consensus_cluster==2~"E2",
                            consensus_cluster==4~"E4",
                            TRUE~NA)) %>% 
  dplyr::rename(a_score=mean_score_01_norm) %>% 
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
peaks_df = peaks %>% 
  filter(distToTSS<=1000) %>% ### do not know if to use the score or the distance from the gene start
  # select(seqnames,start,end,nearestGene,peakType) %>% 
  dplyr::rename(gene=nearestGene) %>% 
  dplyr::left_join(genes_hallmark,relationship = "many-to-many") %>% 
  filter(!is.na(pathway)) %>% 
  mutate(peak=paste(seqnames,start,end,sep='-')) %>% 
  select(peak,pathway) %>% 
  distinct()

pathways = peaks_df$pathway %>% unique()
peak_pathway_list = list()
for (pat in pathways){
  peak_pathway_list[[pat]] <- peaks_df %>% 
    filter(pathway==pat) %>% 
    separate(peak,into = c("chr","from","to"),sep = "-",remove = F,convert = T) %>% 
    mutate(chr = str_remove(chr, "chr"))
}
saveRDS(object = peak_pathway_list,file = 'input_data_P05/data/final_version/peak_pathway_list.rds')
