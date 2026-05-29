library(dplyr)
library(ggplot2)
library(tidyverse)
#library(ComplexHeatmap)
# library(circlize)
library(data.table)
library(parallel)
source("utils.R")
set.seed(12345)
args <- commandArgs(trailingOnly = TRUE)
cell_idx <- as.numeric(args[1])
# cells_interval <- as.numeric(args[1])
# cell_range <- eval(parse(text = cells_interval))

selected_frags_dist_len = readRDS("chr1_selected_frags_dist_len.rds")
df_peak_final = readRDS("df_peak_final_big.rds")
df_peak_cna_final = readRDS("df_peak_cna_final_big.rds")
cell_info = readRDS("cell_info_big.rds")
# selected_cells =df_peak_final$cell_id %>% unique()
# selected_cells =selected_cells[cell_range]
all_cells = df_peak_final$cell_id %>% unique() %>%as.character() %>%  sort()
selected_cells=as.numeric(all_cells[cell_idx])
print(selected_cells)

start=Sys.time()




frag_res_all_cells <- mclapply(
  selected_cells,
  mc.cores = detectCores() - 1,
  function(c){
    
    peak_cell <- df_peak_final %>% 
      filter(cell_id == c) %>% 
      filter(status == 1) %>% 
      separate(
        peak_id,
        into = c("chr","from","to"),
        sep = "-",
        remove = FALSE
      )
    
    peak_cna_cell <- df_peak_cna_final %>% 
      inner_join(peak_cell)
    
    chromosomes <- peak_cna_cell %>% 
      pull(chr) %>% 
      unique()
    
    chromosomes_frags <- lapply(chromosomes, function(chrom){
      
      peak_cna_cell_chr <- peak_cna_cell %>% 
        filter(chr == chrom)
      
      peak_cell_chr <- peak_cell %>% 
        filter(chr == chrom)
      
      fragm_df_cell <- sample_fragments_for_peak_vec(
        peak_id   = peak_cell_chr$peak_id,
        peak_chr=str_replace(chrom,pattern = "chr",replacement = ""),
        peak_from = as.numeric(peak_cell_chr$from),
        peak_to   = as.numeric(peak_cell_chr$to),
        fragment_len_dist = selected_frags_dist_len,
        tot_cn    = as.numeric(peak_cna_cell_chr$tot_cn),
        cell_id   = c
      ) %>% 
        as.data.frame()
      
    }) %>% 
      bind_rows()
    
  }
) %>% 
  bind_rows()
end <- Sys.time()

end - start

### get cell info into peaks

frag_res_all_cells <- frag_res_all_cells %>% inner_join(cell_info)
saveRDS(object = frag_res_all_cells,file = paste0('cell_',selected_cells,'_all_fragments.rds'))
###### create bed file
dir.create(path = "fragments_cells_big")
sampled_cells = cell_info %>% filter(!is.na(sample)) %>% pull(cell_id) %>% unique()

# genome_coordinates = readRDS('genome_coordinates_hg38.rds')
bed_files = lapply(selected_cells,function(c){
  bed_single_cell =frag_res_all_cells %>% 
    filter(cell_id==c) %>% 
    mutate(fragment_start=round(fragment_start,0)) %>% 
    mutate(fragment_end=round(fragment_end,0)) %>% 
    mutate(fragment_id=paste(fragment_chr,fragment_start,fragment_end,allele,peak_id,sep=":")) %>%
    dplyr::select(fragment_chr,fragment_start,fragment_end,fragment_id)
  print(c)
  # chr_start <- genome_coordinates$from[genome_coordinates$chr == "chr22"]
  # bed_single_cell$start_rel <- bed_single_cell$fragment_start - chr_start + 1
  # bed_single_cell$end_rel   <- bed_single_cell$fragment_end   - chr_start + 1
  # bed_single_cell = bed_single_cell %>% dplyr::select(chr,start_rel,end_rel,fragment_id)
  write.table(x = bed_single_cell,file = paste0("fragments_cells_big/cell_",c,".bed"),append = F,quote = F,sep = "\t",row.names = F,col.names = F)
})


########## plot heatmap ###########



# starting_dist_frag_sizes <- replicate(
#   100000,
#   sample_fragment_size()
# )
# 
# library(dplyr)
# library(ggplot2)
# 
# ggplot() +
# 
#   # simulated distribution
#   geom_histogram(
#     data = frag_res_all_cells,
#     aes(x = fragment_size, y = after_stat(density)),
#     binwidth = 1,
#     fill = "steelblue",
#     alpha = 0.5
#   )+
# 
#   # input distribution
#   # geom_histogram(
#   #   data = selected_frags %>% filter(fragment_len <= 800),
#   #   aes(x = fragment_len, y = after_stat(density)),
#   #   binwidth = 1,
#   #   fill = "tomato",
#   #   alpha = 0.5
#   # ) +
#   geom_density(
#     data = data.frame(fragment_len=selected_frags_dist_len[which(selected_frags_dist_len<800)]),
#     aes(x = fragment_len),
#     binwidth = 1,
#     # fill = "tomato",
#     color="tomato",
#     alpha = 0.3
#   )+theme_minimal()
# 

# 
# 
# ####################################
# 
# mutant_cols <- c(
#   "A" = "goldenrod",
#   "B" = "magenta4"
#   # "Clone 3" ="forestgreen"
# )
# 
# epistate_cols <- c(
#   "+" = "forestgreen",
#   "-" = "darkblue"
#   # "Clone 3" ="forestgreen"
# )
# 
# mutant_cols <- c(
#   "A+" = "goldenrod",
#   "A-" = "magenta4",
#   "B+" = "royalblue3",
#   "B-" ="forestgreen"
# )
# peak_class <- c(
#   clonal = "violet",
#   fluctuating = "royalblue3",
#   fixed = "darkorange"
# )
# frag_res_all_cells_all_status <- frag_res_all_cells_all_status %>% 
#   mutate(status=1) %>% 
#   # filter(!str_starts(label.peak_id, "peak_")) %>% 
#   dplyr::select(cell_id,peak_id,mutant,sample,epistate,status) 
# row_annot_df <- frag_res_all_cells_all_status %>%
#   mutate(clone=paste0(mutant,epistate)) %>% 
#   distinct(cell_id, sample,clone) %>%
#   arrange(sample)
# # col_annot_df <- df_final %>% 
# #   distinct(peak_id,type)
# frag_res_all_cells_all_status <- frag_res_all_cells_all_status %>% 
#   mutate(cell_id = factor(cell_id, levels = row_annot_df$cell_id)) %>%
#   arrange(cell_id)
# # convert to matrix: rows = cell_id, cols = peak_id
# mat <- frag_res_all_cells_all_status %>%
#   select(cell_id, peak_id, status) %>%
#   unique() %>% 
#   pivot_wider(
#     names_from = peak_id,
#     values_from = status,values_fill = 0
#   ) %>%
#   column_to_rownames("cell_id") %>%
#   as.matrix()
# 
# # # row annotation dataframe
# # row_annot_mut_df <- final_df %>%
# #   distinct(cell_id, mutant)# %>%
# # #  arrange(cell_id)
# # 
# # row_annot_sample_df <- final_df %>%
# #   distinct(cell_id, sample) #%>%
# #   # arrange(cell_id)
# 
# 
# # colors for mutant annotation
# 
# 
# # row annotation
# 
# 
# row_ha <- rowAnnotation(
#   Clone = row_annot_df$clone,
#   # Epistate = row_annot_df$epistate,
#   # sample =row_annot_df$sample,
#   col = list(Clone=mutant_cols)
# )
# 
# # col_ha <- columnAnnotation(
# #   Peak_Class = col_annot_df$type,
# #   col = list(Peak_Class=peak_class)
# # )
# 
# # heatmap colors
# col_fun <- c(
#   "0" = "white",
#   "1" = "grey"
#   
# )
# 
# # draw heatmap
# 
# 
# plot_heatmap_chromtin <-Heatmap(
#   mat,
#   name = "Chromatin status",
#   col = col_fun,
#   cluster_rows = F,
#   cluster_columns = F,
#   # top_annotation = col_ha,
#   left_annotation = row_ha,show_row_dend = F,
#   show_column_names = F, show_row_names = F#, split = row_annot_df$sample
# )
# # png(filename = "mets_plot/plot_heatmap_chromatin.png",
# #     width = 6, height = 6, units = "in", res = 300)
# draw(plot_heatmap_chromtin)
# # dev.off()
# 
# 
# plot_forest(sample_forest,color_map = mutant_cols)
# 
