suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

cells <- readRDS('cell_info_big.rds')
sampled_cells <- cells %>% 
  filter(!is.na(sample))
fragm_dir <- "fragments_cells_big_with_background_01_lambda_sparsity_070_filtered_peaks/"
rds_files <- list.files(path = fragm_dir,pattern = "rds")
sampled_cells_sequenced <-sapply(str_split(rds_files, "_"), `[`, 2)

simulated_fragments <- lapply(sampled_cells_sequenced, function(cell){
  file_rds <- paste0(fragm_dir,'cell_',cell,'_all_fragments.rds')
  if (file.exists(file_rds)){
    readRDS(file_rds)
  } else {
    message(paste0("RDS for cell ",cell,"does not exists"))
  }
}) %>% bind_rows()

outdir_nextflow<-'analysis_bg/results_v4/'
called_peaks <- read.table(file = paste0(outdir_nextflow,'call_peaks/macs3/sample_A/sample_A_peaks.narrowPeak'), header = F, sep = "\t")
colnames(called_peaks) <- c(
  "chrom", "start", "end", "name", "score",
  "strand", "signalValue", "pValue", "qValue", "peak"
)
called_fragments <- read.table(file=paste0(outdir_nextflow,'fragment_calling/sample_A/sample_A.fragments.tsv.gz'))
colnames(called_fragments) <- c(
  "chrom", "start", "end", "cell", "count"
)
called_fragments=called_fragments %>% 
  mutate(size=end-start) %>% 
  select(!count) %>% 
  mutate(type='inferred')

truth_peaks_bulk_abudance <-simulated_fragments %>% 
  select(peak,cell_id) %>% 
  distinct() %>% 
  group_by(peak) %>% 
  summarise(pct_cells=n()/length(sampled_cells_sequenced))
  
truth_peaks <- simulated_fragments %>% 
  select(peak,peak_chr,peak_from,peak_to) %>% 
  distinct() %>% 
  dplyr::rename(start=peak_from,
         end=peak_to,
         chrom=peak_chr) %>% 
  inner_join(truth_peaks_bulk_abudance)


truth_fragments <- simulated_fragments %>% 
  select(fragment,fragment_chr,fragment_end,fragment_start,fragment_size,cell_id) %>% 
  dplyr::rename(start=fragment_start,
         end=fragment_end,
         chrom=fragment_chr, 
         size=fragment_size) %>% 
  mutate(cell=paste0("cell_",cell_id)) %>% 
  select(chrom,start,end,cell,size) %>% 
  mutate(type='simulated')


fragment_all_long = rbind(truth_fragments,called_fragments)

fragment_size_plot =fragment_all_long %>% 
  ggplot(aes(x=size,fill=type))+
  geom_density(alpha=0.2)+
  theme_minimal()

fragment_count_plot =fragment_all_long %>% 
  ggplot(aes(x=type,fill=type))+geom_bar()+
  theme_minimal()




#-----------------------------------------------------------
# Peak overlap evaluation
#-----------------------------------------------------------

evaluate_peaks <- function(truth, called, min_overlap = 1) {
  truth_gr <- GRanges(
    seqnames = truth$chrom,
    ranges = IRanges(start = truth$start, end = truth$end))
  
  called_gr <- GRanges(
    seqnames = called$chrom,
    ranges = IRanges(start = called$start, end = called$end)
  )
  
  hits <- findOverlaps(
    called_gr,
    truth_gr,
    minoverlap = min_overlap
  )
  colnames(truth)<- paste0('simulated_',colnames(truth))
  colnames(called)<- paste0('inferred_',colnames(called))
  
  matched <- data.frame(
    truth_idx  = subjectHits(hits),
    called_idx = queryHits(hits)
  )
  
  matched <- cbind(
    truth[matched$truth_idx, ],
    called[matched$called_idx, ]
  )

  TP_called <- n_distinct(matched$inferred_name)
  TP_truth  <- n_distinct(matched$simulated_peak_id)
  
  # Peak-level TP count
  TP <- TP_called
  
  # Total peaks
  N_called <- nrow(called)
  N_truth  <- nrow(truth)
  
  
  # False positives and false negatives
  FP <- N_called - TP_called
  FN <- N_truth - TP_truth
  
  # Metrics
  precision <- TP / (TP + FP)
  recall    <- TP / (TP + FN)
  
  f1 <- ifelse(
    precision + recall == 0,
    0,
    2 * precision * recall / (precision + recall)
  )
  
  peak_metrics <- data.frame(
    TP = TP,
    FP = FP,
    FN = FN,
    Precision = precision,
    Recall = recall,
    F1 = f1
  )
  
  return(list('peak_matching_df'=matched,'metrics'=peak_metrics))
  
}


evaluate_peaks_new <- function(truth, called, min_overlap = 1) {
  
  truth_gr <- GRanges(
    seqnames = truth$chrom,
    ranges = IRanges(start = truth$start, end = truth$end)
  )
  
  called_gr <- GRanges(
    seqnames = called$chrom,
    ranges = IRanges(start = called$start, end = called$end)
  )
  
  hits <- findOverlaps(
    called_gr,
    truth_gr,
    minoverlap = min_overlap
  )
  
  colnames(truth)  <- paste0("simulated_", colnames(truth))
  colnames(called) <- paste0("inferred_", colnames(called))
  
  ## -------------------------
  ## True positives (matches)
  ## -------------------------
  matched_idx <- data.frame(
    truth_idx  = subjectHits(hits),
    called_idx = queryHits(hits)
  )
  
  tp_df <- cbind(
    truth[matched_idx$truth_idx, , drop = FALSE],
    called[matched_idx$called_idx, , drop = FALSE]
  )
  tp_df$status <- "TP"
  
  ## -------------------------
  ## False negatives
  ## -------------------------
  fn_idx <- setdiff(seq_len(nrow(truth)), unique(subjectHits(hits)))
  
  fn_df <- truth[fn_idx, , drop = FALSE]
  fn_df$status <- "FN"
  
  ## -------------------------
  ## False positives
  ## -------------------------
  fp_idx <- setdiff(seq_len(nrow(called)), unique(queryHits(hits)))
  
  fp_df <- called[fp_idx, , drop = FALSE]
  fp_df$status <- "FP"
  
  ## Make dataframes compatible for rbind()
  all_cols <- union(colnames(tp_df), union(colnames(fn_df), colnames(fp_df)))
  
  add_missing <- function(df, cols) {
    miss <- setdiff(cols, colnames(df))
    df[miss] <- NA
    df[, cols, drop = FALSE]
  }
  
  tp_df <- add_missing(tp_df, all_cols)
  fn_df <- add_missing(fn_df, all_cols)
  fp_df <- add_missing(fp_df, all_cols)
  
  final_df <- rbind(tp_df, fn_df, fp_df)
  
  ## -------------------------
  ## Metrics
  ## -------------------------
  TP_called <- dplyr::n_distinct(tp_df$inferred_name)
  TP_truth  <- dplyr::n_distinct(tp_df$simulated_peak_id)
  
  TP <- TP_called
  
  N_called <- nrow(called)
  N_truth  <- nrow(truth)
  
  FP <- N_called - TP_called
  FN <- N_truth - TP_truth
  
  precision <- ifelse(TP + FP == 0, 0, TP / (TP + FP))
  recall    <- ifelse(TP + FN == 0, 0, TP / (TP + FN))
  
  f1 <- ifelse(
    precision + recall == 0,
    0,
    2 * precision * recall / (precision + recall)
  )
  
  peak_metrics <- data.frame(
    TP = TP,
    FP = FP,
    FN = FN,
    Precision = precision,
    Recall = recall,
    F1 = f1
  )
  
  return(list(
    peak_matching_df = tp_df,
    peak_classification_df = final_df,
    metrics = peak_metrics
  ))
}

compare_res = evaluate_peaks(truth = truth_peaks,called = called_peaks,min_overlap = 10)
compare_res = evaluate_peaks_new(truth = truth_peaks,called = called_peaks,min_overlap = 10)

all_peaks_classification <- compare_res$peak_classification_df
quantiles <- quantile(compare_res$peak_classification_df$simulated_pct_cells,na.rm=T)[2:4]
all_peaks_classification <- all_peaks_classification %>% 
  mutate(peak_bulk_cell_class=case_when(simulated_pct_cells>=0.8~'high bulk',
                                        simulated_pct_cells<0.8 & simulated_pct_cells>+0.4~'medium bulk',
                                        TRUE~'low bulk'))

all_peaks_classification %>% 
  group_by(peak_bulk_cell_class,status) %>% 
  summarise(n=n())

matching_df=compare_res$peak_matching_df
matching_df <-matching_df %>% 
  mutate(base_difference_start=(simulated_start-inferred_start)) %>% 
  mutate(base_difference_end=(simulated_end-inferred_end)) %>% 
  mutate(abs_diff=abs(base_difference_end)+abs(base_difference_start))

abs_diff_plt=matching_df %>% 
  ggplot(aes(y=abs_diff,x=""))+geom_boxplot(outliers = F)+
  theme_minimal()


truth_peaks_sel=truth_peaks %>% 
  mutate(type="simulated") %>% 
  select(chrom,start,end,type)
called_peaks_sel=called_peaks %>% 
  mutate(type="inferred") %>% 
  select(chrom,start,end,type)
long_df = rbind(called_peaks_sel,truth_peaks_sel)
long_df=long_df %>% 
  mutate(peak_size=end-start)

peak_size_dist = long_df %>% 
  ggplot(aes(y=peak_size,x=type,fill=type))+geom_boxplot(outliers = F)+
  theme_minimal()

peak_counts = long_df %>% 
  ggplot(aes(x=type,fill=type))+geom_bar()+
  theme_minimal()




peak_positions <- data.frame(
  start_diff = matching_df$base_difference_start,
  end_diff   = matching_df$base_difference_end
)

peak_positions <- tidyr::pivot_longer(peak_positions, cols = everything())

peak_position_plt=ggplot(peak_positions, aes(x = value, fill = name)) +
  geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
  labs(
    x = "Difference (bp)",
    y = "Count"
  ) +
  theme_minimal()


performance_metrics = compare_res$metrics
performance_metrics_long <- pivot_longer(
  performance_metrics,
  cols = everything(),
  names_to = "metric",
  values_to = "value"
)
performance_metrics_plt=performance_metrics_long %>% 
  filter(metric%in%c('Precision','Recall','F1')) %>% 
  ggplot(aes(x=metric,y=value, group = 1,fill=metric))+geom_point()+
  geom_col()+
  geom_line()+
  ylim(0,1)+
  theme_minimal()

wrap_plots(list(fragment_count_plot,fragment_size_plot,peak_size_dist,peak_counts,performance_metrics_plt,peak_position_plt,abs_diff_plt),
           guides = 'collect',design = c("AABB\nCCDD\nEFFG"),tag_level = "keep")+plot_annotation(tag_levels = 'a') & theme(legend.position = 'bottom')
