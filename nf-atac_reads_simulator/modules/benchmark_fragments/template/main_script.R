#!/usr/bin/env Rscript
library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(ggplot2)
library(data.table)


source("${projectDir}/bin/utils.R")

# called_fragments <- fread(
#   # "${fragments}",
#   "/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project/analysis/results_case2_genes_5pathways_low_variability/fragment_calling/S3/fragments.chr.tsv.gz",
#   header = FALSE,
#   col.names = c("chrom", "start", "end", "cell_id", "count")
# ) %>%
#   mutate(
#     cell_id = as.character(cell_id),
#     chrom = as.character(chrom),
#     start = as.integer(start),
#     end = as.integer(end)
#   )

output_dir <- "fragment_benchmark"
min_overlap_bp <- 10L
min_overlap_fraction <- 0.50
truth_start_shift <- 0L
truth_end_shift <- 0L
length_breaks <- c(0, 50, 100, 150, 200, 300, 500, Inf)
length_labels <- c("<=50", "51-100", "101-150", "151-200", "201-300", "301-500", ">500")

if (file.exists(file.path(project_dir, "bin/utils.R"))) {
  source(file.path(project_dir, "bin/utils.R"))
}
if (!exists("my_ggplot_theme", mode = "function")) {
  my_ggplot_theme <- function() theme_classic(base_size = 12)
}


# ---- 2. Load called fragments --------------------------------------------
called_file <- "/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project/analysis/results_case2_genes_5pathways_low_variability/fragment_calling/S3/fragments.chr.tsv.gz"
stopifnot(file.exists(called_file))
called_fragments <- fread(called_file, header = FALSE,
                          col.names = c("chrom", "start", "end", "cell_id", "count"))
called <- called_fragments %>%
  transmute(cell_id = as.character(cell_id),
            chrom = as.character(chrom),
            start = as.integer(start), end = as.integer(end),
            count = as.numeric(count)) %>%
  filter(!is.na(cell_id), !is.na(chrom), !is.na(start), !is.na(end),
         end > start, start >= 0)

# ---- 3. Load post-dropout simulated fragments -----------------------------
stopifnot(length(simulated_frags_files) > 0,
          all(file.exists(simulated_frags_files)))
sequenced_fragments_list <- vector("list", length(simulated_frags_files))
for (i in seq_along(simulated_frags_files)) {
  file <- simulated_frags_files[i]
  id <- sub("^fragments_cell_id_([0-9]+)\\.rds$", "\\1", basename(file))
  if (!grepl("^[0-9]+$", id)) stop("Cannot parse cell ID: ", file)
  frags <- readRDS(file)
  sequenced_fragments_list[[i]] <- frags %>%
    filter(sequenced == 1) %>%
    transmute(cell_id = id,
              chrom = as.character(fragment_chr),
              start = as.integer(round(fragment_start)) + truth_start_shift,
              end = as.integer(round(fragment_end)) + truth_end_shift,
              fragment_type = as.character(fragment_type))
}
truth_fragments <- bind_rows(sequenced_fragments_list)
truth <- truth_fragments %>%
  mutate(chrom = if_else(grepl("^chr", chrom), chrom, paste0("chr", chrom))) %>%
  filter(!is.na(cell_id), !is.na(chrom), !is.na(start), !is.na(end),
         end > start, start >= 0)
called <- called %>% mutate(chrom = if_else(grepl("^chr", chrom), chrom, paste0("chr", chrom)))

# Sinto's fifth column is multiplicity; benchmark unique intervals.
# Detect conflicting types before collapsing duplicates.
conflicting_types <- truth %>%
  group_by(cell_id, chrom, start, end) %>%
  summarise(n_types = n_distinct(fragment_type), .groups = "drop") %>%
  filter(n_types > 1)
if (nrow(conflicting_types)) warning("Some identical truth intervals have multiple fragment types; first type retained")
truth <- truth %>% distinct(cell_id, chrom, start, end, .keep_all = TRUE)
called <- called %>% distinct(cell_id, chrom, start, end, .keep_all = TRUE)
cat("Truth unique fragments:", nrow(truth), "\nCalled unique fragments:", nrow(called), "\n")

# ---- 4. Candidate overlaps ------------------------------------------------
# Both datasets assumed BED-like: 0-based, half-open [start, end).
# Preserve global row indices
truth$truth_idx <- seq_len(nrow(truth))
called$called_idx <- seq_len(nrow(called))

# Split by cell
truth_by_cell <- split(truth, truth$cell_id)
called_by_cell <- split(called, called$cell_id)

# Only cells present in both datasets can have matches
common_cells <- intersect(
  names(truth_by_cell),
  names(called_by_cell)
)

matches_list <- vector("list", length(common_cells))

for (i in seq_along(common_cells)) {
  
  cell <- common_cells[i]
  
  truth_cell <- truth_by_cell[[cell]]
  called_cell <- called_by_cell[[cell]]
  
  truth_gr <- GRanges(
    seqnames = truth_cell$chrom,
    ranges = IRanges(
      start = truth_cell$start + 1L,
      end = truth_cell$end
    )
  )
  
  called_gr <- GRanges(
    seqnames = called_cell$chrom,
    ranges = IRanges(
      start = called_cell$start + 1L,
      end = called_cell$end
    )
  )
  
  hits <- findOverlaps(
    truth_gr,
    called_gr,
    minoverlap = min_overlap_bp,
    ignore.strand = TRUE
  )
  
  if (length(hits) == 0L) next
  
  matches_cell <- data.frame(
    truth_idx = truth_cell$truth_idx[queryHits(hits)],
    called_idx = called_cell$called_idx[subjectHits(hits)]
  ) %>%
    mutate(
      cell_id = cell,
      
      # Ground-truth coordinates
      truth_chrom = truth$chrom[truth_idx],
      truth_start = truth$start[truth_idx],
      truth_end = truth$end[truth_idx],
      fragment_type = truth$fragment_type[truth_idx],
      
      # Inferred coordinates
      called_chrom = called$chrom[called_idx],
      called_start = called$start[called_idx],
      called_end = called$end[called_idx],
      
      # Fragment lengths
      truth_length = truth_end - truth_start,
      called_length = called_end - called_start,
      
      # Overlap statistics
      overlap_bp = pmin(truth_end, called_end) -
        pmax(truth_start, called_start),
      
      overlap_truth = overlap_bp / truth_length,
      overlap_called = overlap_bp / called_length,
      
      jaccard = overlap_bp /
        (truth_length + called_length - overlap_bp),
      
      # Coordinate differences
      start_error = called_start - truth_start,
      end_error = called_end - truth_end,
      length_error = called_length - truth_length
    ) %>%
    filter(
      overlap_truth >= min_overlap_fraction,
      overlap_called >= min_overlap_fraction
    )
  matches_list[[i]] <- matches_cell
}
matches <- bind_rows(matches_list) %>%
  arrange(
    desc(jaccard),
    desc(overlap_bp),
    truth_idx,
    called_idx
  )
# ---- 5. Greedy one-to-one matching ----------------------------------------
# Greedy is deterministic, but does not guarantee maximum-cardinality matching.
used_truth <- rep(FALSE, nrow(truth))
used_called <- rep(FALSE, nrow(called))
keep <- logical(nrow(matches))
for (i in seq_len(nrow(matches))) {
  t <- matches$truth_idx[i]
  c <- matches$called_idx[i]
  if (!used_truth[t] && !used_called[c]) {
    keep[i] <- TRUE
    used_truth[t] <- TRUE
    used_called[c] <- TRUE
  }
}





matched_fragments <- matches[keep, , drop = FALSE] %>%
  mutate(start_error = called$start[called_idx] - truth$start[truth_idx],
         end_error = called$end[called_idx] - truth$end[truth_idx],
         length_error = called_length - truth_length)
tp_truth <- truth[matched_fragments$truth_idx, , drop = FALSE]
tp_called <- called[matched_fragments$called_idx, , drop = FALSE]
fn <- truth[!used_truth, , drop = FALSE]
fp <- called[!used_called, , drop = FALSE]

# ---- 6. Overall metrics ---------------------------------------------------
TP <- nrow(matched_fragments); FP <- nrow(fp); FN <- nrow(fn)
precision <- if (TP + FP > 0) TP / (TP + FP) else NA_real_
recall <- if (TP + FN > 0) TP / (TP + FN) else NA_real_
f1 <- if (!is.na(precision) && !is.na(recall) && precision + recall > 0)
  2 * precision * recall / (precision + recall) else NA_real_
metrics <- data.frame(TP = TP, FP = FP, FN = FN,
                      precision = precision, recall = recall, f1 = f1,
                      false_discovery_rate = if (TP + FP > 0) FP / (TP + FP) else NA_real_,
                      mean_jaccard = if (TP > 0) mean(matched_fragments$jaccard) else NA_real_,
                      mean_overlap_truth = if (TP > 0) mean(matched_fragments$overlap_truth) else NA_real_,
                      mean_overlap_called = if (TP > 0) mean(matched_fragments$overlap_called) else NA_real_)
metrics_long <- metrics %>% pivot_longer(everything(), names_to = "metric", values_to = "value")
print(metrics_long)

# ---- 7. Truth recovery annotations ---------------------------------------
truth_eval <- truth %>%
  mutate(truth_idx = row_number(), fragment_length = end - start,
         recovered = used_truth) %>%
  left_join(matched_fragments %>%
              dplyr::select(truth_idx, called_idx, overlap_bp, overlap_truth,
                            truth_chrom,called_chrom,
                            truth_start,called_start,
                            truth_end,called_end,
                     overlap_called, jaccard, start_error, end_error, length_error),
            by = "truth_idx") %>%
  mutate(length_bin = cut(fragment_length, breaks = length_breaks,
                          labels = length_labels, include.lowest = TRUE))

# ---- 8. Recovery by fragment type ----------------------------------------
metrics_by_type <- truth_eval %>%
  group_by(fragment_type) %>%
  summarise(n_truth = n(), TP = sum(recovered), FN = sum(!recovered),
            recall = mean(recovered),
            mean_jaccard = if (sum(recovered) > 0) mean(jaccard[recovered]) else NA_real_,
            median_jaccard = if (sum(recovered) > 0) median(jaccard[recovered]) else NA_real_,
            mean_abs_start_error = if (sum(recovered) > 0) mean(abs(start_error[recovered])) else NA_real_,
            mean_abs_end_error = if (sum(recovered) > 0) mean(abs(end_error[recovered])) else NA_real_,
            .groups = "drop")

# ---- 9. Recovery by simulated length -------------------------------------
metrics_by_length <- truth_eval %>%
  group_by(length_bin, .drop = FALSE) %>%
  summarise(n_truth = n(), TP = sum(recovered), FN = sum(!recovered),
            recall = if (n() > 0) mean(recovered) else NA_real_,
            mean_jaccard = if (sum(recovered) > 0) mean(jaccard[recovered]) else NA_real_,
            .groups = "drop")
metrics_type_length <- truth_eval %>%
  group_by(fragment_type, length_bin, .drop = FALSE) %>%
  summarise(n_truth = n(), TP = sum(recovered),
            recall = if (n() > 0) mean(recovered) else NA_real_, .groups = "drop")

# ---- 10. Per-cell and per-chromosome metrics -------------------------------
# TP truth/called pairs are identical in cell and chromosome by design.
cell_truth <- truth_eval %>% count(cell_id, name = "n_truth")
cell_called <- called %>% count(cell_id, name = "n_called")
cell_tp <- tp_truth %>% count(cell_id, name = "TP")
metrics_by_cell <- full_join(cell_truth, cell_called, by = "cell_id") %>%
  full_join(cell_tp, by = "cell_id") %>%
  mutate(across(c(n_truth, n_called, TP), ~replace_na(.x, 0L)),
         FP = n_called - TP, FN = n_truth - TP,
         precision = if_else(n_called > 0, TP / n_called, NA_real_),
         recall = if_else(n_truth > 0, TP / n_truth, NA_real_),
         f1 = if_else(!is.na(precision) & !is.na(recall) & precision + recall > 0,
                      2 * precision * recall / (precision + recall), NA_real_))
chrom_truth <- truth_eval %>% count(chrom, name = "n_truth")
chrom_called <- called %>% count(chrom, name = "n_called")
chrom_tp <- tp_truth %>% count(chrom, name = "TP")
metrics_by_chrom <- full_join(chrom_truth, chrom_called, by = "chrom") %>%
  full_join(chrom_tp, by = "chrom") %>%
  mutate(across(c(n_truth, n_called, TP), ~replace_na(.x, 0L)),
         FP = n_called - TP, FN = n_truth - TP,
         precision = if_else(n_called > 0, TP / n_called, NA_real_),
         recall = if_else(n_truth > 0, TP / n_truth, NA_real_))

# ---- 11. Exploratory logistic models -------------------------------------
# P-values can be overconfident because fragments are clustered within cells.
model_no_interaction <- NULL
model_recovery <- NULL
if (n_distinct(truth_eval$fragment_type) > 1 &&
    n_distinct(truth_eval$recovered) > 1) {
  model_no_interaction <- glm(recovered ~ log1p(fragment_length) + fragment_type,
                              data = truth_eval, family = binomial())
  model_recovery <- glm(recovered ~ log1p(fragment_length) * fragment_type,
                        data = truth_eval, family = binomial())
  print(summary(model_recovery))
  print(anova(model_no_interaction, model_recovery, test = "Chisq"))
}

# ---- 12. Plots ------------------------------------------------------------
p0 <- metrics_long %>%
  filter(metric %in% c("precision", "recall", "f1", "mean_jaccard",
                       "mean_overlap_truth", "mean_overlap_called")) %>%
  ggplot(aes(x = metric, y = value)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = sprintf("%.3f", value)), vjust = -0.5, size = 3.5) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  coord_cartesian(ylim = c(0, 1.08)) +
  labs(title = "Sinto fragment-calling performance", x = NULL, y = "Metric value") +
  my_ggplot_theme() + theme(axis.text.x = element_text(angle = 35, hjust = 1))

p_type_recall <- ggplot(metrics_by_type,
                        aes(x = fragment_type, y = recall, fill = fragment_type)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.3f", recall)), vjust = -0.4) +
  coord_cartesian(ylim = c(0, 1.08)) +
  scale_fill_manual(values = c("peak"="firebrick4","background"="grey"))+
  labs(title = "Recovery by fragment type", x = "Fragment type", y = "Recall") +
  my_ggplot_theme()

p_type_jaccard <- truth_eval %>% filter(recovered) %>%
  ggplot(aes(x = fragment_type, y = jaccard, fill = fragment_type)) +
  geom_boxplot(outlier.alpha = 0.15, show.legend = FALSE) +
  scale_fill_manual(values = c("peak"="firebrick4","background"="grey"))+
  labs(title = "Overlap quality among recovered fragments", x = "Fragment type",
       y = "Jaccard similarity") + my_ggplot_theme()

p_type_errors <- truth_eval %>% filter(recovered) %>%
  dplyr::select(fragment_type, start_error, end_error) %>%
  pivot_longer(c(start_error, end_error), names_to = "boundary", values_to = "error") %>%
  ggplot(aes(x = fragment_type, y = error, fill = fragment_type)) +
  geom_boxplot(outlier.alpha = 0.15, show.legend = FALSE) +
  scale_fill_manual(values = c("peak"="firebrick4","background"="grey"))+
  facet_wrap(~boundary, scales = "free_y") +
  labs(title = "Coordinate errors by fragment type", x = "Fragment type",
       y = "Called minus simulated (bp)") + my_ggplot_theme()

p_length_recall <- ggplot(metrics_by_length %>% filter(n_truth > 0),
                          aes(x = length_bin, y = recall)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", recall)), vjust = -0.4, size = 3) +
  coord_cartesian(ylim = c(0, 1.08)) +
  labs(title = "Recovery by simulated fragment length", x = "Fragment length (bp)",
       y = "Recall") + my_ggplot_theme()

p_type_length <- ggplot(metrics_type_length %>% filter(n_truth > 0),
                        aes(x = length_bin, y = recall, colour = fragment_type, group = fragment_type)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_color_manual(values = c("peak"="firebrick4","background"="grey"))+
  labs(title = "Recovery by length and fragment type", x = "Fragment length (bp)",
       y = "Recall", colour = "Fragment type") + my_ggplot_theme()

p_length_logistic <- ggplot(truth_eval,
                            aes(x = fragment_length, y = as.integer(recovered), colour = fragment_type)) +
  geom_smooth(method = "glm", method.args = list(family = binomial()),
              formula = y ~ log1p(x), se = TRUE) +
  coord_cartesian(xlim = c(0, 600), ylim = c(0, 1)) +
  scale_color_manual(values = c("peak"="firebrick4","background"="grey"))+
  labs(title = "Probability of fragment recovery", x = "Fragment length (bp)",
       y = "Estimated recovery probability", colour = "Fragment type") +
  my_ggplot_theme()

p_length_distribution <- ggplot(truth_eval,
                                aes(x = fragment_length, fill = fragment_type)) +
  geom_density(alpha = 0.35) + coord_cartesian(xlim = c(0, 600)) +
  scale_fill_manual(values = c("peak"="firebrick4","background"="grey"))+
  labs(title = "Simulated fragment length distributions", x = "Length (bp)",
       y = "Density", fill = "Fragment type") + my_ggplot_theme()

# p_recovered_lengths <- ggplot(truth_eval,
#                               aes(x = fragment_length, fill = recovered)) +
#   geom_density(alpha = 0.35) + facet_wrap(~fragment_type) +
#   coord_cartesian(xlim = c(0, 600)) +
#   labs(title = "Recovered and unrecovered fragment lengths", x = "Length (bp)",
#        y = "Density", fill = "Recovered") + my_ggplot_theme()

p_start_error <- ggplot(matched_fragments, aes(x = start_error)) +
  geom_histogram(bins = 60) +
  labs(title = "Sinto start-coordinate error", x = "Called - truth (bp)", y = "Count") +
  my_ggplot_theme()
p_end_error <- ggplot(matched_fragments, aes(x = end_error)) +
  geom_histogram(bins = 60) +
  labs(title = "Sinto end-coordinate error", x = "Called - truth (bp)", y = "Count") +
  my_ggplot_theme()

plots <- list(overall_metrics = p0, recall_by_type = p_type_recall,
              jaccard_by_type = p_type_jaccard, coordinate_errors_by_type = p_type_errors,
              recall_by_length = p_length_recall, recall_by_type_and_length = p_type_length,
              logistic_recovery = p_length_logistic,
              length_distributions = p_length_distribution,
              # recovered_vs_unrecovered = p_recovered_lengths,
              start_error = p_start_error, end_error = p_end_error)
wrap_plots(plots)

benchmark_metrics <- list(
  parameters = list(min_overlap_bp = min_overlap_bp,
                    min_overlap_fraction = min_overlap_fraction,
                    truth_start_shift = truth_start_shift,
                    truth_end_shift = truth_end_shift),
  overall = metrics,
  overall_long = metrics_long,
  by_fragment_type = metrics_by_type,
  by_length = metrics_by_length,
  by_type_and_length = metrics_type_length,
  by_cell = metrics_by_cell,
  by_chromosome = metrics_by_chrom,
  model = model_results
)

rds_path <- file.path(output_dir, "fragment_benchmark_metrics.rds")
saveRDS(benchmark_metrics, rds_path)

pdf_path <- file.path(output_dir, "fragment_benchmark_report.pdf")
grDevices::pdf(pdf_path, width = 10, height = 7, onefile = TRUE)
tryCatch({
  for (nm in names(plots)) {
    print(plots[[nm]])
  }
}, finally = grDevices::dev.off())

cat("Saved metrics:", normalizePath(rds_path), "\n")
cat("Saved report:", normalizePath(pdf_path), "\n")
