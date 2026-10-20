#!/usr/bin/env Rscript
library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(ggplot2)
# library(ProCESS)
library(tidyr)

source("${projectDir}/bin/utils.R")
called_peaks <- read.table("${peaks}", header = FALSE, sep = "\\t")
colnames(called_peaks) <- c(
  "chrom", "start", "end", "name", "score",
  "strand", "signalValue", "pValue", "qValue", "peak"
)

frags_peak_mapping_files <- strsplit("${cell_peak_txt_files}", " ")[[1]]
simulated_frags_files <- strsplit("${cell_frags_rds_files}", " ")[[1]]


a_score <- readRDS("${activity_list_file}")
activity_df <- convert_activity_list(a_score)
gene_rr <- readRDS("${peak_list_file}")
gene_rr_df <- do.call("rbind",gene_rr)
gene_rr_df_wt_activity <- left_join(gene_rr_df,activity_df,relationship = "many-to-many")

# sample_forest <- load_sample_forest("${sample_forest}")
# tot_cells <-sample_forest\$get_samples_info() %>% 
#   dplyr::filter(name=="${meta.sample_id}") %>% 
#   dplyr::pull(tumour_cells)
# sample_forest_phenoptype_prop <- sample_forest\$get_nodes() %>% 
#   dplyr::filter(sample=="${meta.sample_id}") %>% 
#   dplyr::rename(phenotype=epistate) %>% 
#   dplyr::group_by(phenotype) %>% 
#   dplyr::summarise(prop_cells=n()/tot_cells)
# 
# peak_scores <- gene_rr_df_wt_activity %>% 
#   left_join(sample_forest_phenoptype_prop) %>% 
#   group_by(peak, chr, from, to) %>%
#   summarise(
#     recovery_score = sum(activity * prop_cells),
#     .groups = "drop"
#   )
peak_scores <- gene_rr_df_wt_activity %>% 
  group_by(peak, chr, from, to) %>%
  summarise(
    recovery_score = mean(activity),
    .groups = "drop"
  )

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

background_sequenced_fragments <-list()
for (i in 1:length(simulated_frags_files)){
  sequenced_frags = readRDS(simulated_frags_files[i])
  sequenced_frags = sequenced_frags %>% 
    dplyr::mutate(fragment_start=round(fragment_start,0)) %>% 
    dplyr::mutate(fragment_end=round(fragment_end,0))
  background_sequenced_fragments[[i]] = sequenced_frags %>% 
    dplyr::filter(sequenced==1) %>% 
    dplyr::filter(fragment_type=='background')
  print(i)
}
background_sequenced_fragments_df <- do.call('rbind',background_sequenced_fragments)
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
  ) %>% 
  inner_join(peak_scores)

results <- evaluate_peaks_extended(
  truth = sequenced_peaks_df,
  called = called_peaks,
  min_overlap = 10
)



truth_df <- results\$truth_results %>%
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
  aes(x = fragment_bin, y = recovery,group = 1)
) +
  geom_col(fill = "#4472C4") +
  geom_point(size = 2)+
  geom_line(linewidth = 1)+
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
  aes(x = prevalence_bin, y = recovery,group = 1))+
  geom_point(size = 2)+
  geom_line(linewidth = 1)+
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

####
recovery_by_score <- truth_df %>%
  mutate(
    recovery_score_bin = cut(
      recovery_score,
      breaks = seq(0, 2, by = 0.2),
      include.lowest = TRUE
    )
  ) %>%
  group_by(recovery_score_bin) %>%
  summarise(
    n_regions = n(),
    recovery = mean(recovered),
    .groups = "drop"
  )

p2.0 <- ggplot(
  recovery_by_score,
  aes(x = recovery_score_bin, y = recovery,group = 1)
) +
  geom_col(fill = "darkorange2") +
  geom_point(size = 2)+
  geom_line(linewidth = 1)+
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
called_df <- results\$called_results %>%
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
metrics_long <- results\$metrics %>%
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
    subtitle = paste0("Sample ","${meta.sample_id}"),
    x = NULL,
    y = "Metric value"
  ) +
  my_ggplot_theme()
##### false positive peaks
fp <- results\$fp_peaks %>%
  ungroup()

background <- background_sequenced_fragments_df %>%
  ungroup() %>%
  dplyr::filter(
    sequenced == 1,
    fragment_type == "background"
  ) %>% 
  arrange(
    match(
      sub("^chr", "", as.character(fragment_chr)),
      c(as.character(1:22), "X", "Y", "M", "MT")
    ),
    fragment_start,
    fragment_end
  )

harmonise_chr <- function(x) {
  x <- as.character(x)
  ifelse(grepl("^chr", x), x, paste0("chr", x))
}

# Called peaks: BED coordinates -> GRanges
fp_gr <- GRanges(
  seqnames = harmonise_chr(fp\$chrom),
  ranges = IRanges(
    start = fp\$start + 1L,
    end = fp\$end
  )
)

# Background fragments: assumed 1-based, inclusive
background_gr <- GRanges(
  seqnames = harmonise_chr(background\$fragment_chr),
  ranges = IRanges(
    start = background\$fragment_start,
    end = background\$fragment_end
  )
)
all_chr <- union(
  GenomeInfoDb::seqlevels(fp_gr),
  GenomeInfoDb::seqlevels(background_gr)
)

GenomeInfoDb::seqlevels(fp_gr) <- all_chr
GenomeInfoDb::seqlevels(background_gr) <- all_chr

hits <- findOverlaps(
  query = fp_gr,
  subject = background_gr,
  minoverlap = 10L,
  ignore.strand = TRUE
)

# One row per FP peak, including peaks with no background overlap
fp_background_results <- fp %>%
  mutate(
    n_background_fragments = tabulate(
      S4Vectors::queryHits(hits),
      nbins = nrow(fp)
    ),
    overlaps_background = n_background_fragments > 0L,
    background_overlap = factor(
      if_else(
        overlaps_background,
        "Overlaps background",
        "No background overlap"
      ),
      levels = c("No background overlap", "Overlaps background")
    )
  )

# Summary
fp_background_summary <- fp_background_results %>%
  count(background_overlap, name = "n_peaks", .drop = FALSE) %>%
  mutate(
    fraction = if (sum(n_peaks) > 0L) {
      n_peaks / sum(n_peaks)
    } else {
      NA_real_
    }
  )
p5 <- ggplot(
  fp_background_summary,
  aes(x = background_overlap, y = fraction, fill = background_overlap)
) +
  geom_col(width = 0.65) +
  geom_text(
    aes(
      label = if_else(
        is.na(fraction),
        "No FP peaks",
        paste0(
          scales::percent(fraction, accuracy = 0.1),
          "\n(n = ", n_peaks, ")"
        )
      )
    ),
    vjust = -0.3,
    size = 4
  ) +
  scale_fill_manual(
    values = c(
      "No background overlap" = "#B8BEC6",
      "Overlaps background" = "#3A9D7A"
    ),
    guide = "none"
  ) +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.15))
  ) +
  labs(
    x = NULL,
    y = "Fraction of false-positive peaks",
    title = "Background-fragment overlap with false-positive peaks",
    subtitle = "Overlap requires at least 10 shared base pair"
  ) +
  my_ggplot_theme()

# Include peaks with zero matching fragments
fragment_count_summary <- fp_background_results %>%
  count(n_background_fragments, name = "n_fp_peaks")

p6 <- ggplot(
  fragment_count_summary,
  aes(x = n_background_fragments, y = n_fp_peaks)
) +
  geom_col(width = 0.8, fill = "#3A9D7A") +
  scale_x_continuous(
    breaks = scales::breaks_pretty(),
    minor_breaks = NULL
  ) +
  scale_y_continuous(
    breaks = scales::breaks_pretty(),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = "Number of overlapping background fragments per FP peak",
    y = "Number of false-positive peaks",
    title = "Background-fragment overlap per false-positive peak",
    subtitle = "Overlap requires at least 10 shared base pair"
  ) +
  my_ggplot_theme()



# Store all plots in a named list
final_plots <- list(
  p0 = p0,
  p1 = p1,
  p2 = p2,
  p2.0 = p2.0,
  p3 = p3,
  p4 = p4,
  p5 = p5,
  p6=p6
)


saveRDS(object = results,file = paste0("${meta.sample_id}","_metrics.rds"))
saveRDS(object = final_plots,file = paste0("${meta.sample_id}","_metrics_plots.rds"))

