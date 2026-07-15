
pathway_summary <- final_data %>%
  group_by(consensus_cluster, pathway) %>%
  summarise(
    mean_score = mean(score),
    median_score = median(score),
    sd_score = sd(score),
    n = n(),
    .groups = "drop"
  )


pathway_filter <- pathway_summary %>%
  group_by(pathway) %>%
  summarise(
    max_abs = max(abs(mean_score), na.rm = TRUE),
    sd_across_clusters = sd(mean_score, na.rm = TRUE),
    .groups = "drop"
  )
selected_pathways =pathway_filter %>%
  filter(max_abs > 1) %>%   # adjust threshold (0.7???1.2 typical)
  pull(pathway)

final_data =final_data %>% 
  # filter(pathway%in%selected_pathways) %>% 
  group_by(pathway,consensus_cluster) %>%
  mutate(
    score_01 = scales::rescale(score, to = c(0, 1))
  ) %>%
  ungroup()


pathway_summary %>% 
  ggplot(aes(x=pathway,y=mean_score,color=as.factor(consensus_cluster)))+
  geom_point()+
  # facet_wrap(~class,scales = "free")+
  theme_minimal()+
  theme(axis.text.x = element_blank())


final_data %>% 
  ggplot(aes(x=pathway,y=score,fill=as.factor(consensus_cluster)))+
  geom_boxplot(outliers = F)+
  facet_wrap(~class,scales = "free")+
  theme_minimal()+
  theme(axis.text.x = element_blank())





cluster_mean <- final_data %>%
  group_by(pathway,consensus_cluster) %>%
  summarise(
    mean_score = mean(score),
    sd_score = sd(score)
  )



res <- final_data %>%
  group_by(pathway) %>%
  group_modify(~{
    clusters <- unique(.x$consensus_cluster)
    
    bind_rows(lapply(clusters, function(cl){
      
      in_cl  <- .x$score[.x$consensus_cluster == cl]
      out_cl <- .x$score[.x$consensus_cluster != cl]
      
      tibble(
        cluster = cl,
        mean_in = mean(in_cl),
        mean_out = mean(out_cl),
        effect_size = mean(in_cl) - mean(out_cl)
      )
    }))
  })
top_pathways <- res %>%
  filter(effect_size > 0) %>%
  group_by(cluster) %>%
  arrange(desc(effect_size), .by_group = TRUE) %>%
  slice_head(n = 10)

res %>% 
  left_join(pathway_map_df) %>% 
  # filter(effect_size <= -1.5 | effect_size >= 1.5) %>% 
  # pull(pathway) %>% unique() %>% length()
  ggplot(aes(y=pathway,x=effect_size,color=as.factor(cluster)))+
  geom_point()+
  theme_minimal()
most_variable_pathways = res %>% filter(effect_size <= -1.5 | effect_size >= 1.5) %>% 
  pull(pathway) %>% unique()
mat_effect_pca <- res %>%
  select(pathway, cluster, effect_size) %>%
  pivot_wider(
    names_from = cluster,
    values_from = effect_size
  ) %>%
  column_to_rownames("pathway") %>%
  as.matrix()

mat_effect_pca_filtered =mat_effect_pca[most_variable_pathways,]
pca <- prcomp(mat_effect_pca_filtered, center = TRUE, scale. = TRUE)
pca_res = pca$x %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'cluster')
ggplot(pca_res,aes(x=PC1,y=PC2,color=cluster))+
  geom_point()+
  theme_minimal()+theme(legend.position = "bottom")

mat_effect <- res %>%
  select(pathway, cluster, effect_size) %>%
  mutate(cluster=as.factor(cluster)) %>% 
  mutate(effect_size=as.numeric(effect_size)) %>% 
  tidyr::pivot_wider(
    names_from = cluster,
    values_from = effect_size
  ) %>%
  column_to_rownames("pathway") %>%
  as.matrix()

row_classes <- final_data %>%
  select(pathway, class) %>% 
  distinct()
# tibble::column_to_rownames("pathway")
row_ha <- rowAnnotation(
  class = row_classes$class,
  col = list(
    class = c(
      immune = "#E41A1C",
      stromal = "#377EB8",
      oncogenic = "#4DAF4A",
      stress = "#FF7F00",
      other = "grey70"
    )
  )
)

top_path_ht=Heatmap(mat_effect,left_annotation = row_ha)

draw(
  top_path_ht,
  heatmap_legend_side = "bottom"
)#### Heatmap drawing,

mat <- final_data_filtered %>%
  select(cell_id, pathway, score) %>%
  pivot_wider(names_from = cell_id, values_from = score)
rownames(mat) <- mat$pathway
mat <- mat %>% select(-pathway)
mat <- as.matrix(mat)

row_classes <- final_data_filtered %>%
  select(pathway, class) %>% 
  distinct()
# tibble::column_to_rownames("pathway")
row_ha <- rowAnnotation(
  class = row_classes$class,
  col = list(
    class = c(
      immune = "#E41A1C",
      stromal = "#377EB8",
      oncogenic = "#4DAF4A",
      stress = "#FF7F00",
      other = "grey70"
    )
  )
)

cell_cluster <- final_data_filtered %>%
  select(cell_id, consensus_cluster) %>%
  mutate(consensus_cluster=as.factor(consensus_cluster)) %>% 
  distinct()
# tibble::column_to_rownames("cell_id")
# cell_order <- intersect(colnames(mat), rownames(cell_cluster))

# mat <- mat[, cell_order]
# cluster_vec <- cell_cluster[cell_order, "consensus_cluster"]
column_ha <- HeatmapAnnotation(
  cluster = cell_cluster$consensus_cluster,
  col = list(cluster=c('1'="olivedrab",'2'="olivedrab3",'3'="forestgreen",'4'="palegreen1")
  ))
ht=Heatmap(
  mat,
  name = "Score",
  show_row_names = TRUE,
  show_column_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = T,
  # column_split = cluster_vec,
  
  top_annotation = column_ha,
  left_annotation = row_ha,
  # col = colorRampPalette(c("white", "darkgreen"))(10),
  # col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
)
pdf(file = "heatmap_01_score_filtered.pdf",width = 20,height = 15)
draw(ht)
dev.off()
