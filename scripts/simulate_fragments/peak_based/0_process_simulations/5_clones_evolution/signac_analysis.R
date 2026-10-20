rm(list = ls())

library(Signac)
library(Seurat)
library(GenomicRanges)
library(GenomeInfoDb)
library(Rsamtools)
library(ProCESS)
library(AnnotationHub)
library(dplyr)
library(ggplot2)
library(patchwork)

set.seed(123)

# ------------------------------------------------------------
# 1. Paths and parameters
# ------------------------------------------------------------

project_dir <- paste0(
  "/orfeo/cephfs/scratch/cdslab/ggandolfi/",
  "Github/scATAC_project"
)

source(file.path(
  project_dir,
  "ProCESS-scATAC/scripts/simulate_fragments/peak_based/utils.R"
))

forest_path <- file.path(
  project_dir,
  paste0(
    "ProCESS-scATAC/scripts/simulate_fragments/peak_based/",
    "0_process_simulations/5_clones_evolution/sample_forest.sff"
  )
)

fragment_paths <- setNames(
  file.path(
    project_dir,
    "analysis/results_case5_clones/fragment_calling",
    c("S1", "S2", "S3"),
    "fragments.tsv.gz"
  ),
  c("S1", "S2", "S3")
)

macs_path <- paste0(
  "/orfeo/cephfs/home/cdslab/ggandolfi/.cache/R/",
  "basilisk/1.18.0/zellkonverter/1.16.0/",
  "zellkonverterAnnDataEnv-0.10.9/bin/macs3"
)

simulated_pathway_activities <- readRDS("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project/reference_data/a_scores_gene_case5_clones.rds")
activity_df <-convert_activity_list(simulated_pathway_activities) %>% 
  dplyr::rename(phenotype=epistate) %>% 
  dplyr::rename(gene=pathway)
peaks <- readRDS("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project/reference_data/peak_gene_case5_clones.rds") %>% bind_rows()

fragment_cutoff <- 10
min_features <- 200
min_cells <- 10

stopifnot(
  all(file.exists(fragment_paths)),
  file.exists(forest_path),
  file.exists(macs_path)
)

# ------------------------------------------------------------
# 2. Simulation metadata
# ------------------------------------------------------------

sample_forest <- load_sample_forest(forest_path)

metadata_simulated_cells <- sample_forest$get_nodes() %>%
  dplyr::filter(!is.na(sample)) %>%
  mutate(
    sample = as.character(sample),
    cell_id = as.character(cell_id),
    cell_key = paste(sample, cell_id, sep = "_")
  )

if (anyDuplicated(metadata_simulated_cells$cell_key)) {
  stop("Simulation metadata contains duplicate sample/cell_id pairs.")
}

# ------------------------------------------------------------
# 3. Index fragments, select cells, create fragment objects
# ------------------------------------------------------------

fragment_counts <- list()
barcodes <- list()
fragment_objects <- list()

for (sample_id in names(fragment_paths)) {
  
  message("Preparing fragments for ", sample_id)
  
  fragpath <- fragment_paths[[sample_id]]
  
  # Input must already be coordinate-sorted and BGZF-compressed
  if (!file.exists(paste0(fragpath, ".tbi"))) {
    Rsamtools::indexTabix(
      file = fragpath,
      format = "bed"
    )
  }
  
  fragment_counts[[sample_id]] <- CountFragments(
    fragments = fragpath
  )
  
  barcodes[[sample_id]] <- fragment_counts[[sample_id]] %>%
    dplyr::filter(frequency_count > fragment_cutoff) %>%
    pull(CB) %>%
    as.character()
  
  if (length(barcodes[[sample_id]]) == 0L) {
    stop("No cells passed the fragment cutoff for ", sample_id)
  }
  
  # Names are object cell names; values are barcodes in the file
  cells_map <- setNames(
    barcodes[[sample_id]],
    barcodes[[sample_id]]
  )
  
  fragment_objects[[sample_id]] <- CreateFragmentObject(
    path = fragpath,
    cells = cells_map
  )
}

# ------------------------------------------------------------
# 4. Call peaks separately and build a shared peak set
# ------------------------------------------------------------

called_peaks <- list()

for (sample_id in names(fragment_objects)) {
  
  message("Calling peaks for ", sample_id)
  
  called_peaks[[sample_id]] <- CallPeaks(
    object = fragment_objects[[sample_id]],
    macs2.path = macs_path,
    name = paste0(sample_id, "_MACS3")
  )
}

# Union of the sample peak sets
common_peaks <- GenomicRanges::reduce(
  unlist(
    GRangesList(called_peaks),
    use.names = FALSE
  ),
  ignore.strand = TRUE
)

# Retain standard chromosomes
common_peaks <- keepSeqlevels(
  common_peaks,
  intersect(
    seqlevels(common_peaks),
    standardChromosomes(common_peaks)
  ),
  pruning.mode = "coarse"
)

common_peaks <- sort(common_peaks)

if (length(common_peaks) == 0L) {
  stop("No peaks remain in the shared peak set.")
}

# ------------------------------------------------------------
# 5. Quantify the same peaks in each sample
# ------------------------------------------------------------

sample_objects <- list()

for (sample_id in names(fragment_objects)) {
  
  message("Quantifying shared peaks for ", sample_id)
  
  counts <- FeatureMatrix(
    fragments = fragment_objects[[sample_id]],
    features = common_peaks,
    cells = barcodes[[sample_id]],
    keep_all_features = TRUE
  )
  
  # Keep all features here so every sample has the same peak set.
  # Apply the feature prevalence filter after merging.
  chrom_assay <- CreateChromatinAssay(
    counts = counts,
    sep = c("-", "-"),
    fragments = fragment_objects[[sample_id]],
    min.cells = 0,
    min.features = min_features,
    genome = "hg38"
  )
  
  object <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks",
    project = sample_id
  )
  
  # Match original fragment barcodes to simulation cell IDs
  original_barcodes <- colnames(object)
  
  metadata_keys <- paste(
    sample_id,
    original_barcodes,
    sep = "_"
  )
  
  idx <- match(
    metadata_keys,
    metadata_simulated_cells$cell_key
  )
  
  if (anyNA(idx)) {
    stop(
      "Unmatched simulation metadata for ", sample_id,
      ". Example barcode(s): ",
      paste(
        head(original_barcodes[is.na(idx)]),
        collapse = ", "
      )
    )
  }
  
  # Preserve original IDs before adding sample prefixes
  object$sample <- sample_id
  object$cell_id <- original_barcodes
  
  object$mutant <- as.character(
    metadata_simulated_cells$mutant[idx]
  )
  
  object$phenotype <- sub(
    "^E",
    "P",
    as.character(metadata_simulated_cells$epistate[idx])
  )
  
  # Signac updates the fragment-to-cell mappings when renaming
  object <- RenameCells(
    object = object,
    add.cell.id = sample_id
  )
  
  sample_objects[[sample_id]] <- object
}

# ------------------------------------------------------------
# 6. Merge samples
# ------------------------------------------------------------

pbmc <- merge(
  x = sample_objects[[1]],
  y = unname(sample_objects[-1]),
  project = "CRC_5clones",
  merge.data = FALSE
)

DefaultAssay(pbmc) <- "peaks"

# Keep peaks detected in at least min_cells across all samples
peak_counts <- GetAssayData(
  object = pbmc,
  assay = "peaks",
  layer = "counts"
)

keep_features <- Matrix::rowSums(peak_counts > 0) >= min_cells

pbmc <- subset(
  pbmc,
  features = rownames(peak_counts)[keep_features]
)

# Cells may lose detected features after the joint peak filter
pbmc <- subset(
  pbmc,
  cells = colnames(pbmc)[pbmc$nFeature_peaks >= min_features]
)

table(pbmc$sample)
table(pbmc$sample, pbmc$phenotype)
table(pbmc$sample, pbmc$mutant)

# ------------------------------------------------------------
# 7. Gene annotations
# ------------------------------------------------------------

ah <- AnnotationHub()

ensdb_v98 <- ah[["AH75011"]]

annotations <- GetGRangesFromEnsDb(
  ensdb = ensdb_v98
)

seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg38"

Annotation(pbmc) <- annotations

# ------------------------------------------------------------
# 8. QC metrics
# ------------------------------------------------------------

pbmc <- NucleosomeSignal(object = pbmc)
pbmc <- TSSEnrichment(object = pbmc)

# Retains the AnnotationHub resource used in your original script.
# Confirm this resource corresponds to your intended hg38 blacklist.
blacklist_regions <- ah[["AH107305"]]

seqlevelsStyle(blacklist_regions) <- "UCSC"

pbmc$blacklist_ratio <- FractionCountsInRegion(
  object = pbmc,
  assay = "peaks",
  regions = blacklist_regions
)

pbmc$nucleosome_group <- ifelse(
  pbmc$nucleosome_signal > 4,
  "NS > 4",
  "NS <= 4"
)
samples_colors <- RColorBrewer::brewer.pal(n = 3,name = "Set2")
names(samples_colors) <- sample_forest$get_samples_info()$name


qc_violin <- VlnPlot(
  object = pbmc,
  features = c(
    "nCount_peaks",
    "TSS.enrichment",
    "blacklist_ratio",
    "nucleosome_signal"
  ),
  group.by = "sample",
  pt.size = 0,
  ncol = 4,
  cols = samples_colors
)

qc_violin_pheno <- VlnPlot(
  object = pbmc,
  features = c(
    "nCount_peaks"
  ),
  group.by = "phenotype",
  pt.size = 0,
  ncol = 1,
  cols = phenotype_colors
)
qc_violin_geno <- VlnPlot(
  object = pbmc,
  features = c(
    "nCount_peaks"
  ),
  group.by = "mutant",
  pt.size = 0,
  ncol = 1,
  cols = genotype_colors
)
qc_violin_gp <- wrap_plots(qc_violin_geno,qc_violin_pheno)+
  patchwork::plot_layout(guides = "collect") &
  my_ggplot_theme() &
  theme(legend.position = "bottom")
qc_violin <- qc_violin +
  patchwork::plot_layout(guides = "collect") &
  my_ggplot_theme() &
  theme(legend.position = "bottom")

pbmc@meta.data <- pbmc@meta.data %>% 
  dplyr::mutate(patient="patient1")
Idents(pbmc) <- "patient"
tss_enrichment_plot_scatter <- DensityScatter(
    pbmc,
    x = "nCount_peaks",
    y = "TSS.enrichment",
    log_x = TRUE,
    quantiles = TRUE
  )

simulated_fragm_hist <- FragmentHistogram(
    object = pbmc,
    )+
  my_ggplot_theme()+
  labs(title="Fragment length distribution",fill="Sample")
  # scale_fill_manual(values=samples_colors)
pbmc <- TSSEnrichment(
  object = pbmc,
  fast = FALSE
)

tss_enrichment_plot <- TSSPlot(
  object = pbmc,
  # group.by = "sample"
) +
  my_ggplot_theme()
  # labs(color="Sample")+
  # scale_color_manual(values=samples_colors)
light_grid <- theme(
  panel.grid.minor = element_blank(),
  panel.grid.major = element_line(
    colour = "grey90",
    linewidth = 0.2
  )
)
tss_enrichment_plot< - tss_enrichment_plot+light_grid
simulated_fragm_hist <- simulated_fragm_hist+light_grid
plot_qc <- wrap_plots(
  list(
    wrap_elements(full = qc_violin),
    wrap_elements(full = qc_violin_gp),
    tss_enrichment_plot_scatter,
    simulated_fragm_hist,
    tss_enrichment_plot
  ),
  design = "aaaa\naaaa\nbbcc\nddee"
) +
  # plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  # my_ggplot_theme() &
  theme(legend.position = "bottom")


ggsave(filename = "plot_qc_signac_case_5clones.pdf",plot = plot_qc,width = 12,height = 12)


gene <- "MYC"

tss <- Signac::GetTSSPositions(
  ranges = Annotation(pbmc)
)

gene_tss <- tss[tss$gene_name %in% gene]
stopifnot(length(gene_tss) > 0)

# Highlight ±50 bp around each annotated TSS
highlight <- GenomicRanges::resize(
  gene_tss,
  width = 101,
  fix = "center"
)
highlight$color <- "red"

CoveragePlot(
  object = pbmc,
  region = gene,
  group.by = "sample",
  annotation = TRUE,
  peaks = TRUE,
  region.highlight = highlight,
  extend.upstream = 2000,
  extend.downstream = 2000
)
# ------------------------------------------------------------
# 9. Joint dimensional reduction and clustering
# ------------------------------------------------------------

DefaultAssay(pbmc) <- "peaks"

pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = "q0")
pbmc <- RunSVD(pbmc)

# Use available LSI components, excluding component 1
n_lsi <- ncol(Embeddings(pbmc, reduction = "lsi"))

if (n_lsi < 3L) {
  stop("Too few LSI components for downstream analysis.")
}

lsi_dims <- seq.int(2L, min(30L, n_lsi))

pbmc <- RunUMAP(
  object = pbmc,
  reduction = "lsi",
  dims = lsi_dims
)

pbmc <- FindNeighbors(
  object = pbmc,
  reduction = "lsi",
  dims = lsi_dims
)

pbmc <- FindClusters(
  object = pbmc,
  algorithm = 3,
  verbose = FALSE
)

# ------------------------------------------------------------
# 10. Plot clusters and simulation labels
# ------------------------------------------------------------

phenotype_colors <- c(
  P1 = "forestgreen",
  P2 = "goldenrod",
  P3 = "orchid2"
)

genotype_colors <- c(
  G1 = "coral2",
  G2 = "turquoise4",
  G3 = "darkorange",
  G4 = "mediumorchid3",
  G5 = "limegreen"
)
samples_colors <- RColorBrewer::brewer.pal(n = 3,name = "Set2")
names(samples_colors) <- sample_forest$get_samples_info()$name
umap_clusters <- DimPlot(
  pbmc,
  group.by = "seurat_clusters",
  label = TRUE,
  pt.size = 0.5
) + ggtitle("Signac Clusters")+my_ggplot_theme()

umap_samples <- DimPlot(
  pbmc,
  cols = samples_colors,
  group.by = "sample",
  pt.size = 0.5
) + ggtitle("Samples")+my_ggplot_theme()

umap_phenotypes <- DimPlot(
  pbmc,
  group.by = "phenotype",
  cols = phenotype_colors,
  # label = TRUE,
  pt.size = 0.5
) + ggtitle("Phenotypes")+my_ggplot_theme()

umap_genotypes <- DimPlot(
  pbmc,
  group.by = "mutant",
  cols = genotype_colors,
  # label = TRUE,
  pt.size = 0.5
) + ggtitle("Genetic clones")+my_ggplot_theme()


umap_plots <- wrap_plots(list(umap_clusters,umap_samples,umap_phenotypes,umap_genotypes),
                         design = "aabb\naabb\nccdd\nccdd") + 
  plot_annotation(tag_levels = "A") & theme(legend.position="right") &light_grid
ggsave(filename = "plot_umap_signac_case_5clones.pdf",plot = umap_plots,width = 12,height = 8)


# Replace phenotype and genotype with your metadata column names
md <- pbmc[[]]
table(Cluster = md$seurat_clusters, Phenotype = md$phenotype)
table(Cluster = md$seurat_clusters, Genotype = md$mutant)
table(Cluster = md$seurat_clusters, Sample = md$sample)

mclust::adjustedRandIndex(md$seurat_clusters, md$phenotype)
mclust::adjustedRandIndex(md$seurat_clusters, md$mutant)
mclust::adjustedRandIndex(md$seurat_clusters, md$sample)
### within sample
sapply(split(md, md$sample), function(x) {
  mclust::adjustedRandIndex(x$seurat_clusters, x$phenotype)
})
###à check covrage contribution
lsi <- Embeddings(pbmc, reduction = "lsi")

depth_cor <- apply(lsi, 2, function(x) {
  cor(
    x,
    log1p(pbmc$nCount_peaks),
    method = "spearman"
  )
})

barplot(
  depth_cor,
  las = 2,
  ylab = "Spearman correlation with log1p(nCount_peaks)"
)
# ------------------------------------------------------------
# 11. Gene activity
# ------------------------------------------------------------

DefaultAssay(pbmc) <- "peaks"

gene_activities <- GeneActivity(object = pbmc)

pbmc[["RNA"]] <- CreateAssayObject(
  counts = gene_activities
)

rna_counts <- GetAssayData(
  object = pbmc,
  assay = "RNA",
  layer = "counts"
)

rna_totals <- Matrix::colSums(rna_counts)
positive_totals <- rna_totals[rna_totals > 0]

if (length(positive_totals) == 0L) {
  stop("GeneActivity returned no positive gene-activity counts.")
}

pbmc <- NormalizeData(
  object = pbmc,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = median(positive_totals)
)

DefaultAssay(pbmc) <- "RNA"


pbmc$phenotype <- factor(
  as.character(pbmc$phenotype),
  levels = c("P1", "P2", "P3")
)

pbmc$sample <- factor(pbmc$sample)

stopifnot(
  all(c("RNA", "peaks") %in% Assays(pbmc)),
  !anyNA(pbmc$phenotype),
  !anyNA(pbmc$sample)
)

Idents(pbmc) <- "phenotype"

# Positive effects mean higher activity/accessibility in group 1
contrasts <- list(
  P2_vs_P1 = c("P2", "P1"),
  P3_vs_P1 = c("P3", "P1"),
  P3_vs_P2 = c("P3", "P2")
)

table(pbmc$sample, pbmc$phenotype)
DefaultAssay(pbmc) <- "RNA"

activity_matrix <- GetAssayData(
  pbmc,
  assay = "RNA",
  layer = "data"
)


######## plot heatmap ########
# Log-normalised inferred gene activity
activity_matrix <- GetAssayData(
  pbmc,
  assay = "RNA",
  layer = "data"
)

phenotypes <- c("P1", "P2", "P3")

# Match metadata to matrix column order
cell_phenotypes <- as.character(
  pbmc[[]][colnames(activity_matrix), "phenotype"]
)

# Average across cells within each phenotype
mean_activity <- do.call(
  cbind,
  lapply(phenotypes, function(p) {
    cells <- which(cell_phenotypes == p)
    
    if (length(cells) == 0L) {
      stop("No cells found for ", p)
    }
    
    Matrix::rowMeans(activity_matrix[, cells, drop = FALSE])
  })
)

rownames(mean_activity) <- rownames(activity_matrix)
colnames(mean_activity) <- phenotypes

# Restrict to genes included in the simulation
genes_to_plot <- intersect(
  unique(activity_df$gene),
  rownames(mean_activity)
)

heatmap_matrix <- mean_activity[genes_to_plot, , drop = FALSE]

# Remove genes without variation between phenotypes
gene_sd <- apply(heatmap_matrix, 1, sd)
heatmap_matrix <- heatmap_matrix[
  is.finite(gene_sd) & gene_sd > 0,
  ,
  drop = FALSE
]

# Standardise each gene across phenotypes
heatmap_z <- t(scale(t(heatmap_matrix)))

phenotype_annotation <- rowAnnotation(
  phenotype = colnames(heatmap_z),
  col = list(
    phenotype = phenotype_colors
  )
)

ht <- Heatmap(
  t(heatmap_z),
  name = "Activity\nz-score",
  col = colorRamp2(
    c(-2, 0, 2),
    c("#2166AC", "white", "#B2182B")
  ),
  right_annotation = phenotype_annotation,
  cluster_rows = T,
  cluster_columns = T,
  show_column_names = nrow(heatmap_z) <= 60,
  column_title = "Inferred gene activity by phenotype",
  use_raster = nrow(heatmap_z) > 100
)

draw(ht)

phenotypes <- c("P1", "P2", "P3")

mean_activity <- do.call(
  cbind,
  lapply(phenotypes, function(p) {
    cells <- colnames(pbmc)[pbmc$phenotype == p]
    
    if (length(cells) == 0L) {
      stop("No cells found for phenotype ", p)
    }
    
    Matrix::rowMeans(
      activity_matrix[, cells, drop = FALSE]
    )
  })
)

rownames(mean_activity) <- rownames(activity_matrix)
colnames(mean_activity) <- phenotypes

inferred_activity_long <- as.data.frame(mean_activity) %>%
  rownames_to_column("gene") %>%
  pivot_longer(
    cols = all_of(phenotypes),
    names_to = "phenotype",
    values_to = "inferred_activity"
  ) 



activity_comparison <- inferred_activity_long %>%
  inner_join(
    activity_df,
    by = c("gene", "phenotype")
  ) %>%
  dplyr::filter(
    is.finite(inferred_activity),
    is.finite(activity)
  )

safe_spearman <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  
  if (
    sum(keep) < 3L ||
    length(unique(x[keep])) < 2L ||
    length(unique(y[keep])) < 2L
  ) {
    return(NA_real_)
  }
  
  cor(x[keep], y[keep], method = "spearman")
}

activity_correlations <- activity_comparison %>%
  group_by(phenotype) %>%
  summarise(
    n_genes = n(),
    spearman_rho = safe_spearman(
      activity,
      inferred_activity
    ),
    .groups = "drop"
  )

activity_correlations

p_activity_agreement <- ggplot(
  activity_comparison,
  aes(
    x = activity,
    y = inferred_activity,
    colour = phenotype
  )
) +
  geom_point(alpha = 0.25, size = 0.6) +
  geom_smooth()+
  facet_wrap(~ phenotype) +
  scale_colour_manual(values = phenotype_colors) +
  labs(
    x = "Simulated gene activity",
    y = "Mean inferred gene activity (log-normalised)",
    title = "Agreement between simulated and inferred gene activity"
  ) +
  my_ggplot_theme()+
  theme(legend.position = "none")
p_activity_agreement

############## DAP ##########
DefaultAssay(pbmc) <- "peaks"

pbmc$log10_peak_depth <- log10(pbmc$nCount_peaks + 1)

da_results <- list()

for (comparison in names(contrasts)) {
  
  groups <- contrasts[[comparison]]
  
  message("Testing ", comparison)
  
  object_sub <- subset(
    pbmc,
    cells = colnames(pbmc)[pbmc$phenotype %in% groups]
  )
  
  object_sub$sample <- droplevels(object_sub$sample)
  Idents(object_sub) <- "phenotype"
  
  # Check that phenotype and sample effects can be separated
  design <- model.matrix(
    ~ phenotype + sample,
    data = droplevels(object_sub[[]])
  )
  
  if (qr(design)$rank < ncol(design)) {
    stop(
      comparison,
      ": sample and phenotype are confounded; ",
      "inspect table(pbmc$sample, pbmc$phenotype)."
    )
  }
  
  latent_variables <- "log10_peak_depth"
  
  if (nlevels(object_sub$sample) > 1L) {
    latent_variables <- c(latent_variables, "sample")
  }
  
  markers <- FindMarkers(
    object = object_sub,
    assay = "peaks",
    ident.1 = groups[1],
    ident.2 = groups[2],
    test.use = "LR",
    slot = "counts",
    latent.vars = latent_variables,
    min.pct = 0.05,
    logfc.threshold = 0,
    only.pos = FALSE,
    fc.name = "avg_log2FC",
    verbose = FALSE
  )
  
  da_results[[comparison]] <- markers %>%
    rownames_to_column("peak") %>%
    mutate(
      comparison = comparison,
      significance = case_when(
        p_val_adj < 0.05 & avg_log2FC > 0.25 ~
          paste0("Higher in ", groups[1]),
        
        p_val_adj < 0.05 & avg_log2FC < -0.25 ~
          paste0("Higher in ", groups[2]),
        
        TRUE ~ "Not significant"
      )
    )
}

da_peaks <- bind_rows(da_results)

da_summary <- da_peaks %>%
  count(comparison, significance, name = "n_peaks")

da_summary
########### PLOT DA #############
p_da_volcano <- ggplot(
  da_peaks,
  aes(
    x = avg_log2FC,
    y = -log10(pmax(p_val_adj, .Machine$double.xmin)),
    colour = significance != "Not significant"
  )
) +
  geom_point(alpha = 0.5, size = 0.6) +
  geom_vline(
    xintercept = c(-0.25, 0.25),
    linetype = "dashed",
    colour = "grey60"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    colour = "grey60"
  ) +
  facet_wrap(~ comparison) +
  scale_colour_manual(
    values = c("FALSE" = "grey75", "TRUE" = "#3A9D7A"),
    labels = c("Not significant", "Significant")
  ) +
  labs(
    x = "Log2 fold change in mean peak counts",
    y = expression(-log[10]("Adjusted P-value")),
    colour = NULL,
    title = "Differential accessibility between phenotypes"
  ) +
  theme_classic()

p_da_volcano

significant_peaks <- da_peaks %>%
  dplyr::filter(significance != "Not significant")

if (nrow(significant_peaks) > 0L) {
  
  closest_genes <- ClosestFeature(
    object = pbmc,
    regions = unique(significant_peaks$peak)
  )
  
  da_peaks_annotated <- significant_peaks %>%
    left_join(
      closest_genes,
      by = c("peak" = "query_region")
    )
  
} else {
  
  da_peaks_annotated <- significant_peaks
}
top_peak_p2<-da_peaks_annotated %>% 
  dplyr::filter(significance=="Higher in P3") %>% 
  pull(peak) %>% head(5)
plot2 <- FeaturePlot(
  object = pbmc,
  features = top_peak_p2,
  pt.size = 0.1
)
plot2
