activity <- list(
  'A' = list('+'=c("P1"=0.9, "P2"=0.8, "P3"=0.2),
             '-'=c("P1"=0.1, "P2"=0.3, "P3"=0.3)),
  'B' = list('+'=c("P1"=0.2, "P2"=0.8, "P3"=0.8),
             '-'=c("P1"=0.1, "P2"=0.3, "P3"=0.9))
)

get_epigenetic_activity<- function(activity,mutant,clone){
  programs <- activity[[mutant]][[clone]]
  return(programs)
}

get_epigenetic_activity(activity,mutant="B",clone="-")

# -----------------------------
# 3. Peaks + process assignment
# -----------------------------
n_peaks_P1 <- 100
n_peaks_P2 <- 200
n_peaks_P3 <- 300

peaks <- data.frame(
  peak_id = paste0("peak_", 1:(n_peaks_P1 + n_peaks_P2 + n_peaks_P3)),
  process = c(rep("P1", n_peaks_P1), rep("P2", n_peaks_P2), rep("P3", n_peaks_P3))
)

labelling_functor_new <- function(label, node) {
  
  
  # the nodes are labelled by the identifiers of the associated cells
  cell_id_node = node$cell_id
  cell_mutant = sample_forest$get_nodes() %>%
    filter(cell_id==cell_id_node) %>%
    pull(mutant)
  cell_phenotype = sample_forest$get_nodes() %>%
    filter(cell_id==cell_id_node) %>%
    pull(epistate)
  cell_clone <- paste0(cell_mutant,cell_phenotype)
  clone_programs <-get_epigenetic_activity(activity = activity,mutant = cell_mutant,clone=cell_phenotype)
  cell_activity_peaks <- list()
  for (p in names(clone_programs)){
    program_peaks <- peaks %>% filter(process==p)
    a_score = clone_programs[p]
    
    program_status = rbinom(nrow(program_peaks), size = 1, prob = a_score)
    cell_activity_peaks[[p]]<-peaks %>%
      filter(process==p) %>%
      mutate(status=program_status) %>%
      mutate(cell_id=cell_id_node) %>%
      mutate(mutant=cell_mutant) %>%
      mutate(epistate=cell_phenotype)
  }
  
  final_peaks <- do.call("rbind",cell_activity_peaks)
  
  return(final_peaks)
}
start <- Sys.time()
tour_peaks <- get_label_tour(sample_forest, labelling_functor_new, only_leaves=TRUE)

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