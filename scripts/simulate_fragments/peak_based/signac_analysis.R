rm(list=ls())
library(Signac)
library(Seurat)
library(GenomicRanges)
library(ggplot2)
library(patchwork)
library(GenomeInfoDb)
library(Rsamtools)
library(ProCESS)
library(msigdbr)
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project/ProCESS-scATAC/scripts/simulate_fragments/peak_based/utils.R")
sample_forest <- load_sample_forest("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project/ProCESS-scATAC/scripts/simulate_fragments/peak_based/0_process_simulations/sample_forest_atac_case_1.sff")
metadata_simulated_cells_case1 <- sample_forest$get_nodes() %>% 
  dplyr::filter(!is.na(sample))
fragpath <- '/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project/analysis/results_case1/fragment_calling/S3/fragments.chr.tsv.gz'
indexTabix(
  file = fragpath,
  format = "bed"
)
simulated_pathway_activities <- readRDS("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project/reference_data/a_scores_pathway_case1.rds")
activity_df <-convert_activity_list(simulated_pathway_activities)
peaks <- readRDS("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project/reference_data/peak_pathway_list_unique_peaks.rds")
# Define cells
# If you already have a list of cell barcodes to use you can skip this step
total_counts <- CountFragments(fragpath)
cutoff <- 1000 # Change this number depending on your dataset!
barcodes <- total_counts[total_counts$frequency_count > cutoff, ]$CB

# Create a fragment object
frags <- CreateFragmentObject(path = fragpath, cells = barcodes)

# First call peaks on the dataset
# If you already have a set of peaks you can skip this step
peaks <- CallPeaks(frags,macs2.path = "/orfeo/cephfs/home/cdslab/ggandolfi/.cache/R/basilisk/1.18.0/zellkonverter/1.16.0/zellkonverterAnnDataEnv-0.10.9/bin/macs3")

# Quantify fragments in each peak
counts <- FeatureMatrix(fragments = frags, features = peaks, cells = barcodes)
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = fragpath,
  min.cells = 10,
  min.features = 200
)

pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks"
  # meta.data = metadata
)
peaks.keep <- seqnames(granges(pbmc)) %in% standardChromosomes(granges(pbmc))
pbmc <- pbmc[as.vector(peaks.keep), ]
library(AnnotationHub)
ah <- AnnotationHub()

# Search for the Ensembl 98 EnsDb for Homo sapiens on AnnotationHub
query(ah, "EnsDb.Hsapiens.v98")
ensdb_v98 <- ah[["AH75011"]]
# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = ensdb_v98)

# change to UCSC style since the data was mapped to hg38
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"
# add the gene information to the object
Annotation(pbmc) <- annotations


# compute nucleosome signal score per cell
pbmc <- NucleosomeSignal(object = pbmc)

# compute TSS enrichment score per cell
pbmc <- TSSEnrichment(object = pbmc)

# add fraction of reads in peaks
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100

# add blacklist ratio
blacklist_regions <- ah[['AH107305']] # blacklist regions for hg38
pbmc$blacklist_ratio <- FractionCountsInRegion(
  object = pbmc, 
  assay = 'peaks',
  regions = blacklist_regions
)
DensityScatter(pbmc, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = pbmc, group.by = 'nucleosome_group')


VlnPlot(
  object = pbmc,
  features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'pct_reads_in_peaks'),
  pt.size = 0.1,
  ncol = 5
)
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)
pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)



gene.activities <- GeneActivity(pbmc)
pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)
pbmc <- NormalizeData(
  object = pbmc,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(pbmc$nCount_RNA)
)


# Make sure cell_id is character
metadata_simulated_cells_case1$cell_id <- as.character(metadata_simulated_cells_case1$cell_id)

# Make sure Seurat cell names are character
cell_ids <- rownames(pbmc@meta.data)

# Match the dataframe to the Seurat object
idx <- match(cell_ids, metadata_simulated_cells_case1$cell_id)

# Add metadata columns
pbmc@meta.data <- cbind(
  pbmc@meta.data,
  metadata_simulated_cells_case1[idx, setdiff(colnames(metadata_simulated_cells_case1), "cell_id"), drop = FALSE]
)
DimPlot(object = pbmc, label = TRUE,group.by = "epistate",cols = c("E1"="forestgreen","E2"="goldenrod","E3"="orchid2"),pt.size = 1) 
DimPlot(object = pbmc, label = TRUE,group.by = "mutant",cols = c("G1"="coral2",
                                                                   "G2"="turquoise4",
                                                                   "G3"="darkorange"),pt.size = 1) 

DefaultAssay(pbmc) <- 'RNA'


hallmark <- msigdbr(
  species = "Homo sapiens",
  collection = "H"
)
hallmark <- hallmark %>% 
  dplyr::select(gene_symbol,gs_name) %>% 
  dplyr::rename(gene=gene_symbol) %>% 
  dplyr::rename(pathway=gs_name)

activity_df_with_genes <- activity_df %>% 
  dplyr::full_join(hallmark,relationship = "many-to-many")
activity_df_with_genes_e2 <- activity_df_with_genes %>% 
  dplyr::filter(epistate=="E2") %>% 
  dplyr::filter(pathway=="HALLMARK_KRAS_SIGNALING_DN") %>% 
  dplyr::pull(gene) %>% sample(5)

activity_df_with_genes_e1 <- activity_df_with_genes %>% 
  dplyr::filter(epistate=="E1") %>% 
  dplyr::filter(pathway=="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION") %>% 
  dplyr::pull(gene) %>% sample(5)


activity_df_with_genes_e3 <- activity_df_with_genes %>% 
  dplyr::filter(epistate=="E3") %>% 
  dplyr::filter(pathway=="HALLMARK_DNA_REPAIR") %>% 
  dplyr::pull(gene) %>% sample(10)


# FeaturePlot(
#   object = pbmc,
#   features = c(activity_df_with_genes_e1,activity_df_with_genes_e2,activity_df_with_genes_e3),
#   pt.size = 0.1,
#   max.cutoff = 'q95',
#   ncol=5
# )


FeaturePlot(
  object = pbmc,
  # features = activity_df_with_genes_e3,
  features = c("SMAD5","POLE4","MPG","POLA1","POLD1"),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol=5
)



library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

# Convert to matrix
mat <- activity_df %>%
  dplyr::select(epistate, pathway, activity) %>%
  pivot_wider(
    names_from = epistate,
    values_from = activity
  ) %>%
  column_to_rownames("pathway") %>%
  as.matrix()

# Optional: scale each pathway across epistates

# Heatmap
col_fun <- colorRamp2(c(0, 1), c("white", "steelblue3"))
ht <- Heatmap(
  t(mat),
  name = "Activity",
  
  # Clustering
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  row_dend_reorder = TRUE,
  column_dend_reorder = TRUE,
  
  # Labels
  row_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_gp = gpar(fontsize = 8),
  
  # Cells
  rect_gp = gpar(col = "white", lwd = 0.5),
  
  # Show values
  show_heatmap_legend = TRUE,
  
  # Titles
  column_title = "Pathways",
  col=col_fun
  )

draw(
  ht,
  heatmap_legend_side = "right",
  padding = unit(c(5, 5, 5, 5), "mm")
)



pbmc <- SortIdents(pbmc)
DefaultAssay(pbmc) <- "peaks"
pbmc$epistate <- factor(
  pbmc$epistate,
  levels = c("E1", "E2", "E3")
)

Idents(pbmc) <- "epistate"
p=CoveragePlot(
  object = pbmc,
  region = "SMAD5",
  extend.upstream = 1000,
  extend.downstream = 1000,group.by = "epistate"
)
p &
  scale_fill_manual(
    values = c(
      "E1" = "forestgreen",
      "E2" = "goldenrod",
      "E3" = "orchid2"
    )
  ) 
