library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(ggplot2)

source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project/ProCESS-scATAC/nf-atac_reads_simulator/bin/utils.R")
called_peaks <- read.table("S3_peaks.narrowPeak", header = FALSE, sep = "\t")
colnames(called_peaks) <- c(
  "chrom", "start", "end", "name", "score",
  "strand", "signalValue", "pValue", "qValue", "peak"
)

frags_peak_mapping_files <- list.files(".",pattern = "_peaks_fragment_mapping.txt") 
#strsplit("cell_1812891_peaks_fragment_mapping.txt cell_1287988_peaks_fragment_mapping.txt", " ")[[1]]
simulated_frags_files <- list.files("results_case5_clones/process_simulation/",pattern = "fragments_cell_id_") 
#strsplit("fragments_cell_id_1812891.rds fragments_cell_id_1287988.rds", " ")[[1]]

sequenced_peaks_list<-list()
for (i in 1:length(simulated_frags_files)){
  peak_frags_pre_noise = read.table(frags_peak_mapping_files[i])
  peak_frags_pre_noise = peak_frags_pre_noise %>% 
    dplyr::mutate(fragment_start=round(fragment_start,0)) %>% 
    dplyr::mutate(fragment_end=round(fragment_end,0))
  sequenced_frags = readRDS(simulated_frags_files[i])
  sequenced_frags = sequenced_frags %>% 
    dplyr::mutate(fragment_start=round(fragment_start,0)) %>% 
    dplyr::mutate(fragment_end=round(fragment_end,0))
  sequenced_frags = sequenced_frags %>% 
    filter(sequenced==1) %>% 
    filter(fragment_type=='peak') %>% 
    inner_join(peak_frags_pre_noise,by=c('fragment_start','fragment_end','fragment_allele','fragment_chr','fragment_size'))
  
  sequenced_peaks_list[[i]] = sequenced_frags %>% 
    select(peak, peak_chr,peak_from,peak_to,cell_id)
  print(i)
}

sequenced_peaks_df <- do.call('rbind',sequenced_peaks_list)
tot_cells = length(frags_peak_mapping_files)
sequenced_peaks_df <- sequenced_peaks_df %>%
  ungroup() %>%
  group_by(peak, peak_chr, peak_from, peak_to) %>%
  summarise(
    n_cells = n_distinct(cell_id),
    n_fragments = n(),
    .groups = "drop"
  ) %>%
  mutate(
    pct_cells = n_cells / tot_cells
  )
# sequenced_peaks_df=sequenced_peaks_df %>% 
#   group_by(peak) %>%
#   dplyr::mutate(pct_cells=n()/tot_cells) %>% 
#   dplyr::select(!cell_id) %>% 
#   dplyr::distinct()


# truth=sequenced_peaks_df
# called=called_peaks
# min_overlap=10
# results = evaluate_peaks_new(truth,called,min_overlap)
# saveRDS(object = results,file = paste0("S3","_metrics.rds"))
# 


############# NEW
results <- evaluate_peaks_extended(
  truth = sequenced_peaks_df,
  called = called_peaks,
  min_overlap = 10
)



truth_df <- results$truth_results %>%
  ungroup()
#Recovery as a function of sequenced fragments
recovery_by_fragments <- truth_df %>%
  mutate(
    fragment_bin = cut(
      n_fragments,
      breaks = c(0, 5, 10, 25, 50, 100, 250, 500, Inf),
      labels = c(
        "1-5", "6-10", "11-25", "26-50",
        "51-100", "101-250", "251-500", ">500"
      )
    )
  ) %>%
  group_by(fragment_bin) %>%
  summarise(
    n_regions = n(),
    recovery = mean(recovered),
    .groups = "drop"
  )

p1 <- ggplot(
  recovery_by_fragments,
  aes(x = fragment_bin, y = recovery)
) +
  geom_col(fill = "#4472C4") +
  geom_text(
    aes(label = paste0("n=", n_regions)),
    vjust = -0.3,
    size = 3
  ) +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 1)
  ) +
  labs(
    x = "Sequenced fragments per regulatory region",
    y = "Recovery rate",
    title = "Regulatory-region recovery by fragment abundance"
  ) +
  my_ggplot_theme()
#Recovery as a function of cell prevalence
recovery_by_cells <- truth_df %>%
  mutate(
    prevalence_bin = cut(
      pct_cells,
      breaks = seq(0, 1, by = 0.1),
      include.lowest = TRUE
    )
  ) %>%
  group_by(prevalence_bin) %>%
  summarise(
    n_regions = n(),
    recovery = mean(recovered),
    .groups = "drop"
  )

p2 <- ggplot(
  recovery_by_cells,
  aes(x = prevalence_bin, y = recovery)
) +
  geom_col(fill = "#3A9D7A") +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 1)
  ) +
  labs(
    x = "Fraction of cells with sequenced fragments",
    y = "Recovery rate",
    title = "Regulatory-region recovery by cell prevalence"
  ) +
  my_ggplot_theme()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
# MACS3 peak widths and merging
called_df <- results$called_results %>%
  ungroup() %>%
  mutate(
    peak_width = end - start,
    overlap_category = case_when(
      recovered_truth_regions == 0 ~ "No simulated region",
      recovered_truth_regions == 1 ~ "One simulated region",
      TRUE ~ "Multiple simulated regions"
    )
  )

p3 <- ggplot(
  called_df,
  aes(x = peak_width, fill = overlap_category)
) +
  geom_histogram(
    bins = 60,
    alpha = 0.75,
    position = "identity"
  ) +
  scale_x_log10() +
  facet_wrap(~ overlap_category, scales = "free_y") +
  labs(
    x = "MACS3 peak width (bp, log scale)",
    y = "Number of peaks",
    title = "MACS3 peak width by overlap category"
  ) +
  my_ggplot_theme()+
  theme(legend.position = "none")

#the distribution of simulated regions per called peak

p4 <- called_df %>%
  filter(recovered_truth_regions > 0) %>%
  ggplot(aes(x = recovered_truth_regions)) +
  geom_bar(fill = "#4472C4") +
  scale_x_continuous(breaks = scales::breaks_pretty()) +
  labs(
    x = "Number of simulated regions overlapping a MACS3 peak",
    y = "Number of MACS3 peaks",
    title = "Mapping multiplicity of called peaks"
  ) +
  my_ggplot_theme()

### general stats
metrics_long <- results$metrics %>%
  dplyr::select(precision, recall, f1, jaccard) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = factor(
      Metric,
      levels = c("precision", "recall", "f1", "jaccard"),
      labels = c("Precision", "Recall", "F1-score", "Jaccard")
    )
  )

p0 <- ggplot(metrics_long, aes(x = Metric, y = Value, fill = Metric)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(
    aes(label = sprintf("%.3f", Value)),
    vjust = -0.5,
    size = 4
  ) +
  scale_y_continuous(
    limits = c(0, 1.08),
    breaks = seq(0, 1, 0.2)
  ) +
  labs(
    title = "MACS3 peak-calling performance",
    subtitle = "Sample S3",
    x = NULL,
    y = "Metric value"
  ) +
  my_ggplot_theme()

## Plot an example peak with simulated regions and sequenced fragments
sequenced_fragments_list <- list()

for (i in seq_along(simulated_frags_files)) {
  
  frags <- readRDS(simulated_frags_files[i])
  sequenced_fragments_list[[i]] <- frags %>%
    filter(sequenced == 1) %>%
    transmute(
      cell_id = sub("fragments_cell_id_(\\d+)\\.rds", "\\1", simulated_frags_files[i]),
      fragment_chr = as.character(fragment_chr),
      fragment_start = as.integer(round(fragment_start)),
      fragment_end = as.integer(round(fragment_end)),
      fragment_type
    )
}

sequenced_fragments_df <- bind_rows(
  sequenced_fragments_list
) %>%
  mutate(
    fragment_chr = ifelse(
      grepl("^chr", fragment_chr),
      fragment_chr,
      paste0("chr", fragment_chr)
    )
  )


p5 <- plot_peak_example(
  peak_name = "S3_peak_2",
  results = results,
  fragments = sequenced_fragments_df,
  flank = 300
)
final_plot <- wrap_plots(list(p0,p1,p2,p3,p4),design = "aabbcc\ndddeee")
#############################################################################################
#############################################################################################

