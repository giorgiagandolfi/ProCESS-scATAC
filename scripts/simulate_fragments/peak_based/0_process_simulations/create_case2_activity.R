library(dplyr)
library(tidyr)

set.seed(123)

# your full list of hallmark pathways (example — replace with your actual list)
activity <- readRDS("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project/reference_data/a_scores_pathway.rds")
activity_df <- convert_activity_list(activity_list = activity)
pathways <- activity_df$pathway %>% unique()
hallmark_cCREs <- readRDS("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project/reference_data/cCRE_db/hallmark_pathway_cCRE_25kb.rds")
epistates <- c("E1", "E2", "E3")

###### USE this to test gene level
selected_pathways <- c("HALLMARK_NOTCH_SIGNALING","HALLMARK_MYC_TARGETS_V1",
                       "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION","HALLMARK_WNT_BETA_CATENIN_SIGNALING")

score_matrix <- tribble(
  ~pathway,                                       ~E1,   ~E2,   ~E3,
  "HALLMARK_NOTCH_SIGNALING",                          0.90,  0.15,  0.10,
  "HALLMARK_MYC_TARGETS_V1",                       0.20,  0.90,  0.15,
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",    0.10,  0.15,  0.90,
  "HALLMARK_WNT_BETA_CATENIN_SIGNALING",           0.70,  0.65,  0.70
)

# Reshape to long format: epistate | pathway | score
pathway_scores <- score_matrix %>%
  pivot_longer(cols = -pathway, names_to = "epistate", values_to = "score") %>%
  dplyr::select(epistate, pathway, score) %>%
  arrange(epistate, desc(score))


filtered_cCRE <-hallmark_cCREs %>% 
  dplyr::filter(pathway%in%selected_pathways) %>% 
  dplyr::left_join(pathway_scores, by = "pathway", relationship = "many-to-many")
gene_scores_mean <- filtered_cCRE %>%
  group_by(gene_symbol, epistate) %>%
  summarise(gene_score_pathway = mean(score), .groups = "drop") %>%
  dplyr::mutate(
    gene_score = pmin(pmax(rnorm(n = dplyr::n(), mean = gene_score_pathway, sd = 0.05), 0), 1)
  )
gene_scores_with_ccre <- filtered_cCRE %>%
  dplyr::select(gene_symbol, epistate, ccre_chrom, ccre_start, ccre_end) %>%
  distinct() %>%
  left_join(gene_scores_mean, by = c("gene_symbol", "epistate"))
gene_unique <- gene_scores_with_ccre %>%
  distinct(gene_symbol, epistate, gene_score)



# Build nested named list: list[[epistate]][[gene_symbol]] = score
gene_activity_list <- gene_unique %>%
  split(.$epistate) %>%
  lapply(function(df) {
    setNames(df$gene_score, df$gene_symbol)
  })

saveRDS(object = gene_activity_list,file = "/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project/reference_data/a_scores_gene_case2.rds")
selected_regions <- gene_scores_with_ccre %>% 
  dplyr::select(gene_symbol,ccre_chrom,ccre_start,ccre_end) %>% 
  dplyr::distinct()
genes_unique <- selected_regions$gene_symbol %>% unique()
region_list <- lapply(genes_unique, function(gene){
  selected_regions_gene <- selected_regions %>% 
    dplyr::filter(gene_symbol == gene) %>% 
    dplyr::mutate(
      peak = paste(ccre_chrom, ccre_start, ccre_end, sep = "-")
    ) %>% 
    dplyr::rename(
      pathway = gene_symbol,
      from = ccre_start,
      to = ccre_end
    ) %>% 
    dplyr::mutate(
      chr = gsub("^chr", "", ccre_chrom)
    ) %>% 
    dplyr::select(peak,chr,from,to,pathway)
})
names(region_list)<-genes_unique
saveRDS(object = region_list,file = "/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project/reference_data/peak_gene_list_small.rds")

library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)

# Wide format: genes as rows, epistates as columns
heatmap_matrix <- gene_scores_with_ccre %>%
  distinct(gene_symbol, epistate, gene_score) %>%
  pivot_wider(names_from = epistate, values_from = gene_score) %>%
  tibble::column_to_rownames("gene_symbol") %>%
  as.matrix()

col_fun <- colorRamp2(c(0, 1), c("#F7F7F7","#2166AC"))

Heatmap(
  t(heatmap_matrix),
  name = "Gene score",
  col = col_fun,
  cluster_rows = F,
  cluster_columns = T,
  show_row_names = TRUE,
  show_column_names = F,
  row_names_gp = gpar(fontsize = 8),
  column_title = "Gene activity score across epistates",
  heatmap_legend_param = list(at = c(0, 0.25, 0.5, 0.75, 1))
)



###### USE this when assign randomly
# assign each epistate a distinct, non-overlapping set of "signature" pathways
# n_signature <- sample(8:20, 3, replace = FALSE)
# 
# shuffled <- sample(pathways)
# 
# signature_list <- list(
#   E1 = shuffled[1:n_signature[1]],
#   E2 = shuffled[(n_signature[1] + 1):(n_signature[1] + n_signature[2])],
#   E3 = shuffled[
#     (n_signature[1] + n_signature[2] + 1):
#       (n_signature[1] + n_signature[2] + n_signature[3])
#   ]
# )

# build the dataframe: high activity for signature pathways, low/baseline for the rest
df <- expand_grid(epistate = epistates, pathway = pathways) %>%
  rowwise() %>%
  dplyr::mutate(
    is_signature = pathway %in% signature_list[[epistate]],
    activity = if (is_signature) {
      rnorm(1, mean = 0.85, sd = 0.05)   # high activity for that epistate's pathways
    } else {
      rnorm(1, mean = 0.30, sd = 0.08)   # low/baseline activity otherwise
    }
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(activity = pmin(pmax(activity, 0), 1)) %>%  # clip to [0,1]
  dplyr::select(epistate, pathway, activity)
result <- split(df, df$epistate) %>%
  lapply(function(x) setNames(x$activity, x$pathway))
saveRDS(object = result,file = "/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project/reference_data/a_scores_pathway_case2.rds")
library(tidyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
result <- readRDS("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project/reference_data/a_scores_pathway_case2.rds")
df <- convert_activity_list(result)
wide <- df %>%
  pivot_wider(names_from = epistate, values_from = activity) %>%
  column_to_rownames("pathway")

mat <- as.matrix(wide)

col_fun <- colorRamp2(c(0, 1), c("white", "steelblue3"))
ht=Heatmap(t(mat), col = col_fun, name = "activity",column_names_gp = gpar(fontsize = 8))
pdf("case_2_activity_heatmap.pdf")
draw(ht)
dev.off()
