library(dplyr)
# library(ggplot2)
library(tidyverse)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
# library(data.table)
library(parallel)
source("utils.R")
# set.seed(12345)
args <- commandArgs(trailingOnly = TRUE)
cell_idx <- as.numeric(args[1])
# cells_interval <- as.numeric(args[1])
# cell_range <- eval(parse(text = cells_interval))


########## DESCRIPTION #############
# Input for this script:
# 1. simulated binary matrix
# cell_id                  peak_id status
# 1   54144     chr1-1006219-1006719      1
# 2   54144     chr1-1068981-1069481      1
# 3   54144   chr1-10694533-10695033      1
# 2. single cell tot CN of peak
# cell_id              peak_id tot_cn
# 1   54144 chr1-1006219-1006719      2
# 2   54144 chr1-1068981-1069481      2
# 3   54144 chr1-1136649-1137149      2
# 4   54144 chr1-1157303-1157803      2
# 5   54144 chr1-1158100-1158600      2
# 6   54144 chr1-1162517-1163017      2



selected_frags_dist_len = readRDS("chr1_selected_frags_dist_len.rds")
df_peak_final = readRDS("input_data_P05/df_peak_final_big_new_sparse_085_filtered.rds")
df_peak_cna_final = readRDS("input_data_P05/df_peak_cna_final_big_new_085_filtered.rds")


df_peak_final <- df_peak_final %>% 
  dplyr::select(cell_id,peak,status)
df_peak_cna_final <- df_peak_cna_final %>% 
  dplyr::select(cell_id,peak,tot_cn)


cell_info = readRDS("cell_info_big.rds")
# selected_cells =df_peak_final$cell_id %>% unique()
# selected_cells =selected_cells[cell_range]
all_cells = df_peak_final$cell_id %>% unique() %>%as.character() %>%  sort()
selected_cells=as.numeric(all_cells[cell_idx])
print(selected_cells)

start=Sys.time()


peak_cell <- df_peak_final %>% 
  filter(cell_id == selected_cells) %>% 
  filter(status == 1) %>% 
  separate(
    peak,
    into = c("chr","from","to"),
    sep = "-",
    remove = FALSE
  ) %>% 
  distinct() %>% 
  mutate(from=as.numeric(from),
         to=as.numeric(to))

peak_cna_cell <- df_peak_cna_final %>% 
  filter(cell_id == selected_cells) %>% 
  filter(peak%in%peak_cell$peak) %>% 
  distinct() %>% 
  separate(
    peak,
    into = c("chr","from","to"),
    sep = "-",
    remove = FALSE
  ) %>% 
  mutate(from=as.numeric(from),
         to=as.numeric(to))


chromosomes <- peak_cell %>% 
  pull(chr) %>% 
  unique()


chromosomes_frags <- mclapply(chromosomes, mc.cores = detectCores() - 1,
                              function(chrom){
  message(paste0("chromosome: ",chrom))  
  peak_cna_cell_chr <- peak_cna_cell %>% 
    filter(chr == chrom)
  
  peak_cell_chr <- peak_cell %>% 
    filter(chr == chrom)
  
  fragm_df_cell <- sample_fragments_for_peak_vec(
    peak_id   = peak_cell_chr$peak,
    peak_chr=str_replace(chrom,pattern = "chr",replacement = ""),
    peak_from = as.numeric(peak_cell_chr$from),
    peak_to   = as.numeric(peak_cell_chr$to),
    fragment_len_dist = selected_frags_dist_len,
    tot_cn    = as.numeric(peak_cna_cell_chr$tot_cn),
    cell_id   = selected_cells
  ) %>%  as.data.frame()
}) %>% bind_rows()


end <- Sys.time()

end - start
message("Peak fragments generated")
### extract regions that are not peak
cell_bg_regions <- get_background_regions(peak_fragments_df = chromosomes_frags,genome = 'hg38',
                                          gaps_file = 'gap.txt.gz',filter_small_than = 150,
                                          centromeres_file = 'cytoBand.txt.gz')
final_mapping <- readRDS('genome_representation.rds')
frag_len_out_peak = final_mapping %>% 
  filter(region_type=='out peak') %>% 
  # ggplot(aes(x=fragment_len))+geom_density()
  pull(fragment_len)
frag_len_out_peak_dens <- density(frag_len_out_peak,from=100)


background_frg=simulate_background_fragments(background_regions = cell_bg_regions,
                                             lambda_per_kb = 0.1,frag_len_out_peak_dens = frag_len_out_peak_dens)


message("Background fragments generated")
### get cell info into peaks
outidr <- "fragments_cells_big_with_background_01_lambda_sparsity_085_filtered_peaks_tss/"

dir.create(path = outidr)
chromosomes_frags <- chromosomes_frags %>% inner_join(cell_info)
saveRDS(object = chromosomes_frags,file = paste0(outidr,'cell_',selected_cells,'_all_fragments.rds'))
###### create bed file

sampled_cells = cell_info %>% filter(!is.na(sample)) %>% pull(cell_id) %>% unique()





bed_files = lapply(selected_cells, function(c){

  bed_single_cell_peak <- chromosomes_frags %>% 
    filter(cell_id == c) %>% 
    mutate(
      fragment_start = formatC(round(fragment_start, 0), format = "f", digits = 0),
      fragment_end   = formatC(round(fragment_end, 0), format = "f", digits = 0),
      fragment_id = paste(fragment_chr, fragment_start, fragment_end, allele, peak_id, sep=":")
    ) %>%
    select(fragment_chr, fragment_start, fragment_end, fragment_id)
  bed_single_cell_bg <- background_frg %>% 
    mutate(
      fragment_chr = str_remove(frag_chr,'chr'),
      fragment_start = formatC(round(frag_start, 0), format = "f", digits = 0),
      fragment_end   = formatC(round(frag_end, 0), format = "f", digits = 0),
      fragment_id = paste(fragment_chr, fragment_start, fragment_end,sep=":")
    ) %>%
    select(fragment_chr, fragment_start, fragment_end, fragment_id)
  # convert to data.table (fast)
  bed_single_cell <- rbind(bed_single_cell_peak,bed_single_cell_bg)
  setDT(bed_single_cell)

  fwrite(
    bed_single_cell,
    file = paste0(outidr,"cell_", c, ".bed"),
    sep = "\t",
    col.names = FALSE,
    quote = FALSE,
    nThread = parallel::detectCores()
  )

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
