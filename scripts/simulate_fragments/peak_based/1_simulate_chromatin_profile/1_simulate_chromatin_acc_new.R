# rm(list=ls())
library(ProCESS)
library(dplyr)
library(tidyverse)
library(readxl)
library(parallel)
source("utils.R")

args <- commandArgs(trailingOnly = TRUE)
#### acticity and gene files
activity_list_file <- (args[1])
peak_df_file <- (args[2])


#### phylo forest path
phylo_forest_file <- args[3]


#### outfiles
out_cna_file <- (args[4])
out_peak_acc_file <- (args[5])
##### The main inputs required for this script are:
##### for each of the simulated process an accessibility score
##### for all epigenetic clones
##### activity <- list(
# 'A' = list('+'=c("P1"=0.9, "P2"=0.8, "P3"=0.2),
#            '-'=c("P1"=0.1, "P2"=0.3, "P3"=0.3)),
# 'B' = list('+'=c("P1"=0.2, "P2"=0.8, "P3"=0.8),
#            '-'=c("P1"=0.1, "P2"=0.3, "P3"=0.9))
# )
# dataframe with peaks and associated process
# peak                     pathway   
# <chr>                    <chr>     
#   1 chr1-109763665-109764165 CRC_TISSUE
# 2 chr1-1136649-1137149     CRC_TISSUE
# 3 chr1-118960498-118960998 CRC_TISSUE
# 4 chr1-167082621-167083121 CRC_TISSUE
# 5 chr1-202085778-202086278 CRC_TISSUE
# 6 chr1-20250565-20251065   CRC_TISSUE
# sample_forest <- load_sample_forest("sample_forest_atac_epigenome.sff")
phylo_forest <- load_phylogenetic_forest(phylo_forest_file)

# activity <- list(
#   'A' = list('+'=c("P1"=0.9, "P2"=0.8, "P3"=0.2),
#              '-'=c("P1"=0.1, "P2"=0.3, "P3"=0.3)),
#   'B' = list('+'=c("P1"=0.2, "P2"=0.8, "P3"=0.8),
#              '-'=c("P1"=0.1, "P2"=0.3, "P3"=0.9))
# )



# -----------------------------
# 3. Peaks + process assignment
# -----------------------------
# n_peaks_P1 <- 100
# n_peaks_P2 <- 200
# n_peaks_P3 <- 300
# 
# peaks <- data.frame(
#   peak_id = paste0("peak_", 1:(n_peaks_P1 + n_peaks_P2 + n_peaks_P3)),
#   process = c(rep("P1", n_peaks_P1), rep("P2", n_peaks_P2), rep("P3", n_peaks_P3))
# )

activity <- readRDS(activity_list_file)
peaks <-  readRDS(peak_df_file) %>% 
  dplyr::rename(process=pathway)
  # dplyr::rename(process=gene)

cna_data_sim <- phylo_forest$get_cell_allelic_fragmentation()
cna_data_sim = cna_data_sim %>% 
  mutate(chr=paste0('chr',chr))

labelling_functor_new <- function(label, node) {
  
  
  # the nodes are labelled by the identifiers of the associated cells
  cell_id_node = node$cell_id
  cell_mutant = phylo_forest$get_nodes() %>%
    filter(cell_id==cell_id_node) %>%
    pull(mutant)
  cell_phenotype = phylo_forest$get_nodes() %>%
    filter(cell_id==cell_id_node) %>%
    pull(epistate)
  cell_clone <- paste0(cell_mutant,cell_phenotype)
  clone_programs <-get_epigenetic_activity(activity = activity,mutant = cell_mutant,clone=cell_phenotype)
  print(cell_mutant)
  print(cell_phenotype)
  cell_activity_peaks <- mclapply(names(clone_programs), function(p) {
    
    program_peaks <- peaks %>%
      filter(process == p)
    
    a_score <- clone_programs[[p]]
    
    program_status <- rbinom(
      n = nrow(program_peaks),
      size = 1,
      prob = a_score
    )
    
    program_peaks %>%
      mutate(
        status = program_status,
        cell_id = cell_id_node,
        mutant = cell_mutant,
        epistate = cell_phenotype
      )
    
  }, mc.cores = detectCores() - 1)
  
  names(cell_activity_peaks) <- names(clone_programs)
  
  final_peaks <- bind_rows(cell_activity_peaks, .id = "process_name")
  
  return(final_peaks)
}

labelling_functor3 <- function(label, node) {
  
  
  # the nodes are labelled by the identifiers of the associated cells
  cell_id_node = node$cell_id
  cell_mutant = phylo_forest$get_nodes() %>% 
    filter(cell_id==cell_id_node) %>% 
    pull(mutant)
  cell_phenotype = phylo_forest$get_nodes() %>% 
    filter(cell_id==cell_id_node) %>% 
    pull(epistate)
  
  cell_cna = cna_data_sim %>% filter(cell_id==cell_id_node) %>% 
    dplyr::rename(from=begin) %>% 
    dplyr::rename(to=end) %>% 
    select(chr,major,minor,from,to)
  
  if (nrow(cell_cna)!=0){
    ### update the status dataframe for each peak
    df_peak_types_cell <- peaks %>%
      mutate(cell_id=cell_id_node) %>% 
      separate(peak,into = c('chr','from','to'),sep = '-',remove = F) %>% 
      mutate(from=as.numeric(from),
             to=as.numeric(to)) 
    peaks_dt <- as.data.table(df_peak_types_cell)
    cna_dt   <- as.data.table(cell_cna)
    setkey(peaks_dt, chr, from, to)
    setkey(cna_dt, chr, from, to)
    res =foverlaps(
      peaks_dt,
      cna_dt,
      by.x = c("chr", "from", "to"),
      by.y = c("chr", "from", "to"),
      type = "any",   # overlap allowed
      nomatch = NA
    )
    df_peak_with_cna <- as.data.frame(res) %>% 
      mutate(tot_cn=major+minor) %>% 
      dplyr::select(cell_id,peak,tot_cn)
  } else {
    df_peak_with_cna <- peaks %>%
      mutate(cell_id=cell_id_node) %>% 
      mutate(tot_cn=NA)
  }
  return(df_peak_with_cna)
}

start <- Sys.time()
tour_peaks <- get_label_tour(phylo_forest, labelling_functor_new, only_leaves=TRUE)
tour_cna <- get_label_tour(phylo_forest, labelling_functor3, only_leaves=TRUE)
end <- Sys.time()
end - start

list_peaks_final <- list()
i=1
start <- Sys.time()
while (!tour_peaks$done) {
  # print(tour$value)
  list_peaks_final[[i]] <- tour_peaks$value
  i=i+1
  print(i)
  tour_peaks$step()
}
end <- Sys.time()
end - start
df_peak_final <- do.call("rbind",list_peaks_final)

list_peak_cna_final <- list()
i=1
start <- Sys.time()
while (!tour_cna$done) {
  # print(tour$value)
  list_peak_cna_final[[i]] <- tour_cna$value
  i=i+1
  tour_cna$step()
}
end <- Sys.time()
end - start
df_peak_cna_final <- do.call("rbind",list_peak_cna_final)



df_peak_final <- df_peak_final %>% 
  # dplyr::select(cell_id,label.peak,label.status) %>% 
  dplyr::rename(peak=label.peak) %>% 
  dplyr::rename(status=label.status) #%>% 
  # dplyr::select(cell_id,peak,status)
df_peak_cna_final <- df_peak_cna_final %>% 
  # dplyr::select(cell_id,label.peak,label.tot_cn) %>%
  dplyr::rename(peak=label.peak) %>%
  dplyr::rename(tot_cn=label.tot_cn) #%>%
  # dplyr::select(cell_id,peak,tot_cn)
cell_info <- phylo_forest$get_nodes()

# df_peak_final_dropout <- add_sparsity(real_df = df_peak_final,dropout_rate = 0.85)
# saveRDS(object = df_peak_final_dropout,file = "input_data_P05/df_peak_final_big_new_sparse_085_filtered.rds")

saveRDS(object = df_peak_cna_final,file = out_cna_file)
saveRDS(object = df_peak_final,file = out_peak_acc_file)
print("Done")
# saveRDS(object = cell_info,file = "cell_info_big.rds")
