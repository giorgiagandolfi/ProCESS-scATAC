library(ComplexHeatmap)
library(Matrix)

df_peak_final = readRDS("input_data_P05/data/df_peak_final_big_new_085_filtered.rds")
df_peak_final_dropout = readRDS("input_data_P05/data/df_peak_final_big_new_sparse_085_filtered.rds")
# df_peak_cna_final = readRDS("input_data_P05/df_peak_cna_final_big_new.rds")

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



mat_all <- df_peak_final %>% 
  # filter(label.process%in%c('HALLMARK_APOPTOSIS','HALLMARK_MYC_TARGETS_V1')) %>% 
  mutate(peak_id_complete=paste0(peak,"_",label.process)) %>% 
  select(cell_id,peak_id_complete,status) %>% 
  mutate(status=as.numeric(status)) %>% 
  pivot_wider(
    names_from = peak_id_complete,
    values_from = status,
    values_fill = 0
  )
mat_all <- as.matrix(mat_all[, -1])
rownames(mat_all) <- df_peak_final |> distinct(cell_id) |> pull(cell_id)

mat_all_sparse <- df_peak_final_dropout %>% 
  # filter(label.process%in%c('HALLMARK_APOPTOSIS','HALLMARK_MYC_TARGETS_V1')) %>% 
  mutate(peak_id_complete=paste0(peak,"_",label.process)) %>% 
  select(cell_id,peak_id_complete,status) %>% 
  mutate(status=as.numeric(status)) %>% 
  pivot_wider(
    names_from = peak_id_complete,
    values_from = status,
    values_fill = 0
  )
mat_all_sparse <- as.matrix(mat_all_sparse[, -1])
rownames(mat_all_sparse) <- df_peak_final_dropout |> distinct(cell_id) |> pull(cell_id)


sampled_peaks <- sample(x = colnames(mat_all),size = 5000,replace = F)
mat = mat_all[,sampled_peaks]
mat_sparse = mat_all_sparse[,sampled_peaks]


mutant_cols <- c(
  "A+" = "goldenrod",
  "A-" = "magenta4",
  "B+" = "royalblue3",
  "B-" ="forestgreen"
)


pathway_class_cols <- c(
  'other'='grey',
  'colon-specific' = 'hotpink2',
  'stress'='khaki4',
  'immune'='firebrick4',
  'oncogenic'='blue4',
  'stromal'='forestgreen'
  
)
row_ann_df <- df_peak_final %>% 
  # filter(label.process%in%c('HALLMARK_APOPTOSIS','HALLMARK_MYC_TARGETS_V1')) %>% 
  mutate(peak_id_complete=paste0(peak,"_",label.process)) %>% 
  mutate(geno_epi=paste0(label.mutant,label.epistate)) %>% 
  select(cell_id,geno_epi) %>% 
  distinct()

row_ann_df <- row_ann_df[match(rownames(mat), row_ann_df$cell_id), ]
row_ann <- rowAnnotation(epistate=row_ann_df$geno_epi,col=list(epistate=mutant_cols))

col_ann_df <- df_peak_final %>% 
  # filter(label.process%in%c('HALLMARK_APOPTOSIS','HALLMARK_MYC_TARGETS_V1')) %>% 
  mutate(peak_id_complete=paste0(peak,"_",label.process)) %>% 
  mutate(geno_epi=paste0(label.mutant,label.epistate)) %>% 
  select(peak_id_complete,label.process) %>% 
  dplyr::rename(pathway=label.process) %>% 
  distinct() %>% 
  left_join(pathway_map_df) %>% 
  mutate(
    class = case_when(
      is.na(class) & pathway == "CRC_TISSUE" ~ "colon-specific",
      is.na(class) ~ "other",
      TRUE ~ class
    ))

col_ann_df <- col_ann_df[match(colnames(mat), col_ann_df$peak_id_complete), ]
col_ann <- HeatmapAnnotation(pathway_class=col_ann_df$class,col=list(pathway_class=pathway_class_cols))


ht=Heatmap(
  mat,
  show_row_dend = F,
  cluster_columns = T,
  show_column_dend = F,
  show_column_names = F,show_row_names = F,
  name = "status",
  col = c("0" = "white", "1" = "grey"),
  left_annotation = row_ann,top_annotation = col_ann
)
pdf(file = 'plots/simulated_peaks_no_sparsity.pdf')
draw(ht)
dev.off()

ht=Heatmap(
  mat_sparse,
  show_row_dend = F,
  cluster_columns = T,
  show_column_dend = F,
  show_column_names = F,show_row_names = F,
  name = "status",
  col = c("0" = "white", "1" = "grey"),
  left_annotation = row_ann,top_annotation = col_ann
)
pdf(file = 'simulated_peaks_with_sparsity_095.pdf')
draw(ht)
dev.off()


activity_list <- readRDS('input_data_P05/data/activity_list.rds')
activity_df <- convert_activity_list(activity_list) %>% 
  mutate(epigenetic_class = paste0(mutant, epistate)) %>% 
  left_join(pathway_map_df) %>% 
  mutate(
    class = case_when(
      is.na(class) & pathway == "CRC_TISSUE" ~ "colon-specific",
      is.na(class) ~ "other",
      TRUE ~ class
    )
  )

activity_matrix <- activity_df %>%
  select(epigenetic_class, pathway, activity) %>%
  
  pivot_wider(
    names_from = pathway,
    values_from = activity
  )
epigenetic_classes = (activity_matrix$epigenetic_class)
activity_matrix <- as.matrix(activity_matrix[, -1])
rownames(activity_matrix) <-epigenetic_classes

row_ann <- rowAnnotation(epigenetic_class=rownames(activity_matrix),col=list(epigenetic_class=mutant_cols))
col_ann_df <- activity_df %>% 
  select(pathway,class) %>% 
  distinct()
col_ann <- HeatmapAnnotation(pathway_class=col_ann_df$class,col=list(pathway_class=pathway_class_cols))
ht=Heatmap(activity_matrix,
        cluster_columns = T,
        top_annotation = col_ann,
        right_annotation = row_ann,
        column_names_gp = gpar(fontsize = 8),
        col = colorRampPalette(c("white", "darkgreen"))(10),)
pdf(file = 'plots/simulated_activity_matrix.pdf',width = 15,height = 10)
draw(ht)
dev.off()





activity_list <- readRDS('input_data_P05/data/activity_list_gene_level_log1pscore.rds')
activity_df <- convert_activity_list(activity_list) %>% 
  mutate(epigenetic_class = paste0(mutant, epistate))
  # left_join(pathway_map_df) %>% 
  # mutate(
  #   class = case_when(
  #     is.na(class) & pathway == "CRC_TISSUE" ~ "colon-specific",
  #     is.na(class) ~ "other",
  #     TRUE ~ class
  #   )
  # )

activity_matrix <- activity_df %>%
  select(epigenetic_class, pathway, activity) %>%
  
  pivot_wider(
    names_from = pathway,
    values_from = activity
  )
epigenetic_classes = (activity_matrix$epigenetic_class)
activity_matrix <- as.matrix(activity_matrix[, -1])
rownames(activity_matrix) <-epigenetic_classes

row_ann <- rowAnnotation(epigenetic_class=rownames(activity_matrix),col=list(epigenetic_class=mutant_cols))

ht=Heatmap(activity_matrix,
           cluster_columns = T,
           right_annotation = row_ann,
           column_names_gp = gpar(fontsize = 8),
           col = colorRampPalette(c("white", "darkgreen"))(10),)
pdf(file = 'plots/simulated_activity_matrix_per_gene.pdf',width = 15,height = 10)
draw(ht)
dev.off()
