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
sample_forest <- load_sample_forest("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project/ProCESS-scATAC/scripts/simulate_fragments/peak_based/0_process_simulations/5_clones_evolution/sample_forest.sff")
metadata_simulated_cells <- sample_forest$get_nodes() %>% 
  dplyr::filter(!is.na(sample))
fragpath <- '/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project/analysis/results_case5_clones/fragment_calling/S2/fragments.tsv.gz'
indexTabix(
  file = fragpath,
  format = "bed"
)
simulated_pathway_activities <- readRDS("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project/reference_data/a_scores_gene_case5_clones.rds")
activity_df <-convert_activity_list(simulated_pathway_activities)
peaks <- readRDS("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project/reference_data/peak_gene_case5_clones.rds") %>% bind_rows()

# Define cells
# If you already have a list of cell barcodes to use you can skip this step
total_counts <- CountFragments(fragpath)
cutoff <- 10 # Change this number depending on your dataset!
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
# pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100

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
umap_0 <-DimPlot(object = pbmc, label = TRUE,pt.size = 1)

# Make sure cell_id is character
metadata_simulated_cells$cell_id <- as.character(metadata_simulated_cells$cell_id)

# Make sure Seurat cell names are character
cell_ids <- rownames(pbmc@meta.data)

# Match the dataframe to the Seurat object
idx <- match(cell_ids, metadata_simulated_cells$cell_id)

# Add metadata columns
pbmc@meta.data <- cbind(
  pbmc@meta.data,
  metadata_simulated_cells[idx, setdiff(colnames(metadata_simulated_cells), "cell_id"), drop = FALSE]
)
pbmc@meta.data <- pbmc@meta.data %>% 
  dplyr::rename(phenotype=epistate) %>%
  dplyr::mutate(phenotype = stringr::str_replace(phenotype, "^E", "P"))


umap_1 <-DimPlot(object = pbmc, label = TRUE,group.by = "phenotype",cols = c("P1"="forestgreen","P2"="goldenrod","P3"="orchid2"),pt.size = 1) 
genotype_colors <- c("G1"="coral2","G2"="turquoise4", "G3"="darkorange","G4"="mediumorchid3","G5"="limegreen")
umap_2 <-DimPlot(object = pbmc, label = TRUE,group.by = "mutant",cols = genotype_colors,pt.size = 1) 
umap_0+umap_1+umap_2


gene.activities <- GeneActivity(pbmc)
pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)
pbmc <- NormalizeData(
  object = pbmc,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(pbmc$nCount_RNA)
)
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
activity_df_with_genes_p2 <- activity_df %>% 
  dplyr::filter(epistate=="P2") %>% 
  dplyr::slice_max(order_by = activity, n = 5) %>% 
  # dplyr::filter(pathway=="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION") %>% 
  dplyr::pull(pathway) %>% head(5)

activity_df_with_genes_p1 <- activity_df %>% 
  dplyr::filter(epistate=="P1") %>% 
  dplyr::slice_max(order_by = activity, n = 5) %>% 
  # dplyr::filter(pathway=="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION") %>% 
  dplyr::pull(pathway) %>% head(5)


activity_df_with_genes_p3 <- activity_df %>% 
  dplyr::filter(epistate=="P3") %>% 
  # dplyr::filter(pathway=="HALLMARK_DNA_REPAIR") %>% 
  dplyr::slice_max(order_by = activity, n = 5) %>% 
  # dplyr::filter(pathway=="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION") %>% 
  dplyr::pull(pathway) %>% head(5)

ft_plot_ge1 <-FeaturePlot(
  object = pbmc,
  features = c(activity_df_with_genes_p1),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol=5,cols = c("gainsboro","forestgreen")
)
ft_plot_ge2 <-FeaturePlot(
  object = pbmc,
  features = c(activity_df_with_genes_p2),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol=5,cols = c("gainsboro","goldenrod")
)

ft_plot_ge3 <-FeaturePlot(
  object = pbmc,
  features = c(activity_df_with_genes_p3),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol=5,cols = c("gainsboro","orchid2")
)



ft_plot_all <-FeaturePlot(
  object = pbmc,
  features = c("PSEN2","MYC","DVL2","WNT6","JAG2"),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol=5
)

# ft_plot =FeaturePlot(
#   object = pbmc,
#   # features = activity_df_with_genes_e3,
#   features = c("SMAD5","MPG","COL4A1","MYC"),
#   pt.size = 0.1,
#   max.cutoff = 'q95',
#   ncol=4
# )
report_umap <-wrap_plots(list(umap_0,umap_1,umap_2,ft_plot_ge1,ft_plot_ge2,ft_plot_ge3,ft_plot_all),design = "AABBCC\nDDDDDD\nEEEEEE\nFFFFFF\nGGGGGG")
ggsave(filename = "umap_case2_genes_5pathways_low_var.pdf",plot = report_umap,width = 15,height = 15)
pathways_ass_peak=readRDS("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project/reference_data/peak_pathway_list_unique_peaks.rds") %>% 
  bind_rows()


library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

# Convert to matrix
mat <- activity_df %>%
  dplyr::select(epistate, pathway, activity) %>%
  dplyr::rename(phenotype=epistate) %>%
  dplyr::mutate(phenotype = stringr::str_replace(phenotype, "^E", "P")) %>% 
  pivot_wider(
    names_from = phenotype,
    values_from = activity
  ) %>%
  column_to_rownames("pathway") %>%
  as.matrix()
mat <- mat[int_genes,]
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
  show_column_names = F,
  # Labels
  row_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_gp = gpar(fontsize = 4),
  
  # Cells
  # rect_gp = gpar(col = "white", lwd = 0.5),
  
  # Show values
  show_heatmap_legend = TRUE,
  
  # Titles
  column_title = "Genes",
  col=col_fun
  )
pdf("gene_activity_case5clone_int_genes.pdf",width = 20,height = 8)
draw(
  ht,
  heatmap_legend_side = "right",
  padding = unit(c(5, 5, 5, 5), "mm")
)
dev.off()


pbmc <- SortIdents(pbmc)
DefaultAssay(pbmc) <- "peaks"
pbmc$phenotype <- factor(
  pbmc$phenotype,
  levels = c("P1", "P2", "P3")
)

Idents(pbmc) <- "phenotype"
p=CoveragePlot(
  object = pbmc,
  region = "MYC",
  extend.upstream = 1000,
  extend.downstream = 1000,group.by = "phenotype",annotation = TRUE,region.highlight = myc_tss
)
p &
  scale_fill_manual(
    values = c(
      "P1" = "forestgreen",
      "P2" = "goldenrod",
      "P3" = "orchid2"
    )
  ) 





DefaultAssay(pbmc) <- 'peaks'

# wilcox is the default option for test.use
da_peaks <- FindMarkers(
  object = pbmc,
  ident.1 = "0",
  ident.2 = "1",
  test.use = 'wilcox',
  min.pct = 0.1
)
plot1 <- VlnPlot(
  object = pbmc,
  features = rownames(da_peaks)[1:5],
  pt.size = 0.1,
  idents = c("0","1")
)
plot2 <- FeaturePlot(
  object = pbmc,
  features = rownames(da_peaks)[1:5],
  pt.size = 0.1
)

open_c0 <- rownames(da_peaks[da_peaks$avg_log2FC > 3, ])
open_c1 <- rownames(da_peaks[da_peaks$avg_log2FC < -3, ])

closest_genes_c0 <- ClosestFeature(pbmc, regions = open_c0)
closest_genes_c1 <- ClosestFeature(pbmc, regions = open_c1)
