# rm(list=ls())
library(dplyr)
library(purrr)
library(tidyverse)
library(readxl)
library(GenomicRanges)
library(ComplexHeatmap)
library(circlize)
source("utils.R")


cell_epigenetic_df = read.table("/data/rds/DMP/UCEC/GENEVOD/ccolson/MECCA/AllCells/HallmarkPathwayAnalysis/HIPEC_P05/C4_Normalised/magic_gsva_bootstrap_zscoreParam_consensus_clusters.csv",header = T,sep=",")

peaks = readRDS("/data/rds/DMP/UCEC/GENEVOD/ccolson/MECCA/AllCells/ArchRProjects/HIPEC_P05_AllCells/PeakCalls/Type/Tumour-reproduciblePeaks.gr.rds")%>% as.data.frame()
peaks_df = peaks %>% 
  filter(distToTSS<=1000) %>% ### do not know if to use the score or the distance from the gene start
  select(seqnames,start,end,nearestGene,peakType) %>% 
  dplyr::rename(gene=nearestGene)

peaks_per_genes_df = peaks_df %>% 
  dplyr::group_by(gene) %>% 
  dplyr::summarise(tot_peaks=n())

normalised_gene_score <- read.csv("/data/rds/DMP/UCEC/GENEVOD/ccolson/MECCA/AllCells/HallmarkPathwayAnalysis/HIPEC_P05/AllTumourCells_Gene_Normalised_Zscores_10-03-26.csv",row.names = 1)

mat_long <- normalised_gene_score %>% 
  rownames_to_column(var = "cell_id")  %>% 
  pivot_longer(
    cols = -cell_id,
    names_to = "gene",
    values_to = "score"
  ) %>% 
  mutate(score=as.numeric(score)) %>% 
  dplyr::filter(gene%in%peaks_df$gene)
  
  # dplyr::left_join(peaks_per_genes_df) %>% 
  # dplyr::filter(!is.na(tot_peaks)) %>% 
  # dplyr::mutate(normalized_score=score/tot_peaks)

# mat_long %>% 
#   dplyr::group_by(gene) %>% 
#   dplyr::mutate(mean_score=mean(score)) %>% 
#   dplyr::select(gene,mean_score,tot_peaks) %>% 
#   distinct() %>% 
#   ggplot(aes(x=tot_peaks,y=mean_score))+
#   geom_point()

mean_scores <- mat_long %>% 
  dplyr::left_join(cell_epigenetic_df,relationship = "many-to-many") %>%
  dplyr::filter(!is.na(consensus_cluster)) %>%
  dplyr::group_by(gene,consensus_cluster) %>%
  dplyr::summarise(mean_score=mean(score)) %>% 
                   # mean_score_peak_norm=mean(normalized_score)) %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(mean_score_norm =(mean_score - min(mean_score)) /
                  (max(mean_score) - min(mean_score))) %>% 
  # dplyr::mutate(mean_score_norm_per_peak =(mean_score_peak_norm - min(mean_score_peak_norm)) /
  #                 (max(mean_score_peak_norm) - min(mean_score_peak_norm))) %>% 
  dplyr::group_by(gene) %>%
  dplyr::mutate(
    mean_score_norm_per_epi = (mean_score - min(mean_score)) /
      (max(mean_score) - min(mean_score))
  )



# mat_mean_scores_norm_group_long <- mat_mean_scores_norm_group %>% 
#   as.data.frame() %>%
#   rownames_to_column("pathway") %>%
#   pivot_longer(
#     cols = -pathway,
#     names_to = "class",
#     values_to = "score"
#   )
# 
# pathway_peak_df <- peaks_df %>% 
#   dplyr::filter(pathway%in%mat_mean_scores_norm_group_long$pathway)

sampled_genes <- mean_scores %>% pull(gene) %>% unique() %>% sample(50)
sampled_genes =c(sampled_genes,"MYC","CD44","IL6")
test = mean_scores %>% filter(gene%in%sampled_genes)
# test = mean_scores %>% 
#   filter(gene%in%peaks_df$pathway)
# test =test %>% 
#   group_by(gene) %>%
#   mutate(
#     log_score = log1p(mean_score),
#     gene_effect = median(log_score),
#     epi_effect = log_score - gene_effect
#   ) %>%
#   ungroup()
# test$p_open <- plogis(
#   2 * scale(test$gene_effect) +
#     2 * test$epi_effect
# )


test = test %>% 
  group_by(gene) %>%
  mutate(
    score_z = as.numeric(scale(mean_score))
  ) %>%
  ungroup()

test <- test %>%
  group_by(gene) %>%
  mutate(
    gene_baseline = mean(mean_score)
  ) %>%
  ungroup()



test <-test %>%
  mutate(log_baseline = log1p(gene_baseline)) %>% 
  mutate(
    gene_z = as.numeric(scale(log_baseline))
  ) %>% 
  mutate(
    p_open = plogis(
      1.5 * gene_z +
        1.0 * score_z
    )
  )


# test_mat =test %>% 
#   select(gene,consensus_cluster,p_open) %>% 
#   distinct() %>% 
#   tidyr::pivot_wider(
#     names_from = consensus_cluster,
#     values_from = p_open,values_fill = 0
#   ) %>%
#   as.data.frame()
test_mat_raw =test %>% 
  select(gene,consensus_cluster,mean_score) %>% 
  distinct() %>% 
  tidyr::pivot_wider(
    names_from = consensus_cluster,
    values_from = mean_score,values_fill = 0
  ) %>%
  as.data.frame()

test_mat_zscore =test %>% 
  select(gene,consensus_cluster,mean_score_norm) %>% 
  distinct() %>% 
  tidyr::pivot_wider(
    names_from = consensus_cluster,
    values_from = mean_score_norm,values_fill = 0
  ) %>%
  as.data.frame()

rownames(test_mat) <- test_mat$gene
test_mat$gene <- NULL
test_mat <- as.matrix(test_mat)


rownames(test_mat_raw) <- test_mat_raw$gene
test_mat_raw$gene <- NULL
test_mat_raw <- as.matrix(test_mat_raw)


rownames(test_mat_zscore) <- test_mat_zscore$gene
test_mat_zscore$gene <- NULL
test_mat_zscore <- as.matrix(test_mat_zscore)

pdf("normalized_gene_heatmap_new.pdf",height = 20)
ht=Heatmap(
  test_mat,
  name = "p_open",
  cluster_rows = TRUE,
  cluster_columns = F,
  show_row_names = T,
  show_column_names = T,
  # row_names_gp = gpar(fontsize = 3),
  col = colorRampPalette(c("white", "darkgreen"))(10)
)
# draw(ht)


hraw=Heatmap(
  test_mat_raw,
  name = "mean_gene_score",
  cluster_rows = TRUE,
  cluster_columns = F,
  show_row_names = T,
  show_column_names = T,
  # row_names_gp = gpar(fontsize = 3),
  col = colorRampPalette(c("white", "orange"))(20)
)
# draw(hraw)

hz=Heatmap(
  test_mat_zscore,
  name = "mean_zscore",
  cluster_rows = TRUE,
  cluster_columns = F,
  show_row_names = T,
  show_column_names = T,
  col = colorRampPalette(c("white", "darkgreen"))(20)
  # row_names_gp = gpar(fontsize = 3),
  # col = colorRampPalette(c("blue", "white","red"))(10)
)
draw(hz)

ht_list =  hraw + hz
draw(ht_list)
# dev.off()
