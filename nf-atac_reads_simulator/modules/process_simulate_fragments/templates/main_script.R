#!/usr/bin/env Rscript

# Template R script for in-peak fragment simulation
# Nextflow template variables are injected at runtime

# --- Inputs from Nextflow process ---
cell_id<- "${meta.id}"
sample_id<- "${meta.sample_id}"
sample_forest_file<- "${sample_forest}"
phylo_forest_file<- "${pyhlo_forest}"
activity_list_file<- "${activity_list_file}"
peak_list_file<- "${peak_list_file}"
fragment_size_distribution<- "${fragment_size_distribution}"
blacklist_file<- "${encode_blacklist}" #"ENCFF356LFX.bed.gz"
gap_file <- "${gap_file}"
centromere_file <- "${centromere_file}"
fragm_len_dist_out <- "${fragm_len_dist_out}"




library(ProCESS)
library(dplyr)
library(tidyverse)
library(parallel)
library(cli)

source("${projectDir}/bin/utils.R")



sample_forest <- load_sample_forest(sample_forest_file)
phylo_forest <- load_phylogenetic_forest(phylo_forest_file)
reference_genome_path <- phylo_forest\$get_reference_path()
reference_dir <- dirname(reference_genome_path)


activity <- readRDS(activity_list_file)
peaks <-  readRDS(peak_list_file)
selected_frags_dist_len = readRDS(fragment_size_distribution)


peak_access_labelling <- function(node) {
  cell_id_node = node\$cell_id
  message(paste0("Cell: ",cell_id_node))
  
  cell_mutant <- node\$mutant_name
  cell_epistate <- node\$epistate_name
  clone_programs <-get_epigenetic_activity(activity = activity,epistate=cell_epistate)
  clone_programs=sort(clone_programs)
  cell_activity_peaks <- list()
  for (p in names(clone_programs)){
    
    program_peaks <- peaks[[p]]
    a_score = clone_programs[p]
    
    program_status = rbinom(nrow(program_peaks), size = 1, prob = a_score)
    
    cell_activity_peaks[[p]]<-program_peaks %>%
      mutate(
        status = program_status,
        cell_id = cell_id_node,
        mutant = cell_mutant,
        epistate = cell_epistate
      )
    message(paste0("Processing pathway: ",p))
    
  }
  cell_active_peaks <- do.call("rbind",cell_activity_peaks)
  return(cell_active_peaks)
}

node_tour <- get_node_tour(phylo_forest, only_leaves = T,with_genomes = T)

node <- phylo_forest\$get_node(as.numeric(cell_id))
cell_peaks <- peak_access_labelling(node)
saveRDS(object = cell_peaks,file = paste0('cell_',cell_id,'_peak_accessibility.rds'))
genome <- node\$get_genome()
simulated_frags_list <- list()


fasta_cell <- file(paste0("cell_",cell_id,".fasta"), open = "w")

# cell_peaks <- add_sparsity_per_cell(cell_peaks,weibull_scale = 0.85)

for (cp in 1:nrow(cell_peaks)){
  if (cell_peaks\$status[cp]==1){
    cell_alleles = genome\$get_alleles_covering_ref_region(cell_peaks\$chr[cp],
                                                           cell_peaks\$from[cp],
                                                           cell_peaks\$peak_lenght[cp])
    simulated_frags_cell <- sample_fragments_for_peak_vec_allele(
      peak_id   = cell_peaks\$peak[cp],
      peak_chr = cell_peaks\$chr[cp],
      peak_from = as.numeric(cell_peaks\$from[cp]),
      peak_to   = as.numeric(cell_peaks\$to[cp]),
      fragment_len_dist = selected_frags_dist_len,
      available_alleles = cell_alleles,
      cell_id   = cell_id
    ) %>%  as.data.frame() %>%
      ungroup()
    simulated_frags_list[[cp]]<-simulated_frags_cell
    
  } else {
    message("Closed chromatin peak")
  }
}
simulated_frags_df = do.call("rbind",simulated_frags_list)
write.table(x = simulated_frags_df,file = paste0("cell_",cell_id,"_peaks_fragment_mapping.txt",append = F,quote = F,sep = '\t'))


chrom_sizes_file = read.table(file.path(reference_dir,'reference.fasta.chi'))
chrom_sizes_file = chrom_sizes_file %>% dplyr::select(V1,V3)


cell_bg_regions <- get_background_regions(peak_fragments_df = simulated_frags_df,
                                              gaps_file = gap_file,filter_small_than = 150,
                                              centromeres_file = centromere_file,chrom_sizes_file = chrom_sizes_file,
                                              blacklist_file=blacklist_file)


final_mapping <- readRDS(fragm_len_dist_out)
frag_len_out_peak = final_mapping %>% 
  filter(region_type=='out peak') %>% 
  pull(fragment_len)
#frag_len_out_peak_dens <- density(frag_len_out_peak,from=100)
background_frg=simulate_background_fragments(background_regions = cell_bg_regions,
                                             mean_pct_fragments_out = 0.4,
                                             tot_fragments_in_peak = nrow(simulated_frags_df),
                                             frag_len_out_peak_dens = frag_len_out_peak)
background_frg = background_frg %>%
  dplyr::mutate(fragment_type='background') %>% 
  dplyr::mutate(fragment_allele=0)
simulated_frags_df = simulated_frags_df %>% 
  dplyr::mutate(fragment_type='peak') %>%
  select(colnames(background_frg)) 

all_fragments = rbind(simulated_frags_df,background_frg)
all_fragments = all_fragments %>%
  dplyr::mutate(fragment_allele=as.numeric(fragment_allele),
                fragment_start=as.numeric(fragment_start),
                fragment_size=as.numeric(fragment_size))

p_dropout = 0.9
all_fragments\$sequenced <- rbinom(
  n = nrow(all_fragments),
  size = 1,
  prob = 1-p_dropout)
saveRDS(object = all_fragments,file = paste0('fragments_cell_id_',cell_id,".rds"))

n_fragments_to_sequence <- all_fragments %>% 
  dplyr::filter(sequenced==1) %>% 
  nrow()
fragments_to_seq <- all_fragments %>% 
  dplyr::filter(sequenced==1)

pb <- cli_progress_bar(
  "Processing fragments",
  total = n_fragments_to_sequence,
  format = "Processing fragments {cli::pb_bar} {cli::pb_percent} ({cli::pb_current}/{cli::pb_total})"
)
for (f in seq_len(n_fragments_to_sequence)){
  if (fragments_to_seq\$fragment_type[f]=='background'){
    cell_alleles_backgound = genome\$get_alleles_covering_ref_region(fragments_to_seq\$fragment_chr[f],
                                                                     fragments_to_seq\$fragment_start[f],
                                                                     fragments_to_seq\$fragment_size[f])
    random_selected_allele <- sample(x = cell_alleles_backgound,size = 1)
    fragment_region <- genome\$get_region_aligned_to_ref(fragments_to_seq\$fragment_chr[f],
                                                         random_selected_allele %>% as.numeric(),
                                                         fragments_to_seq\$fragment_start[f],
                                                         fragments_to_seq\$fragment_size[f])
    fragment = genome\$get_fragment(fragment_region\$chr,
                                    fragment_region\$allele,
                                    fragment_region\$from,
                                    fragment_region\$length)
  } else {
    fragment_region <- genome\$get_region_aligned_to_ref(fragments_to_seq\$fragment_chr[f],
                                                         fragments_to_seq\$fragment_allele[f],
                                                         fragments_to_seq\$fragment_start[f],
                                                         fragments_to_seq\$fragment_size[f])
    
    fragment = genome\$get_fragment(fragment_region\$chr,
                                    fragment_region\$allele,
                                    fragment_region\$from,
                                    fragment_region\$length)
  }

  fragment_seq = fragment\$sequence
  fragment_id = paste0(">",cell_id,":",fragment_region\$chr,":", fragment_region\$from,":",fragment_region\$allele,":",fragments_to_seq\$fragment_type[f])
  writeLines(fragment_id, fasta_cell)
  writeLines(fragment_seq, fasta_cell)
  cli_progress_update()
}
cli_progress_done()
