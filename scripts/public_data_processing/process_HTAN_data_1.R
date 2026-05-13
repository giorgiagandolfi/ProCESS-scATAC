############################################
## scATAC PRIMARY vs METASTASIS ANALYSIS ##
## Signac full pipeline                   ##
############################################

library(Seurat)
library(Signac)
library(Rsamtools)
library(GenomicRanges)
library(GenomeInfoDb)
library(harmony)

# library(JASPAR2020)
# library(TFBSTools)
# library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)
# setwd("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC/")

############################################
## 1. INPUT DATA
############################################

# Metastasis tumour
barcodes_met <- readLines("CM618C2-T1-barcodes.tsv")
peaks_met <- read.table("CM618C2-T1-peaks.bed", sep = "\t")
counts_met  <- Matrix::readMM("CM618C2-T1-matrix.mtx")
peaks_met_gr <- GRanges(
  seqnames = peaks_met$V1,
  ranges = IRanges(start = peaks_met$V2, end = peaks_met$V3)
)
peaknames_met <- paste0(peaks_met$V1, ":",peaks_met$V2,"-",peaks_met$V3)
colnames(counts_met) <- barcodes_met
rownames(counts_met) <- peaknames_met


fragpath_met <- 'CM618C2-T1-fragments.tsv.gz'

rds_met <- readRDS("CRC_CM618C1-T1Y2.rds")
metadata_met <- rds_met@meta.data
  
chrom_assay_met <- CreateChromatinAssay(
  counts = counts_met,
  sep = c(":", "-"),
  fragments = fragpath_met,
  min.cells = 10,
  min.features = 200
)

met_obj <- CreateSeuratObject(
  counts = chrom_assay_met,
  assay = "peaks",
  meta.data = metadata_met
)

met_obj <- subset(
  met_obj,
  cells = rds_met$Original_barcode
)
cell_type_vector_names <- gsub(pattern = "CRC_CM618C1-T1Y2_",replacement = "",x = names(rds_met$cell_type))
cell_type_vector <- rds_met$cell_type
names(cell_type_vector) <- cell_type_vector_names
met_obj$cell_type <- cell_type_vector





# Primary tumour

barcodes_primary <- readLines("CM618C1-S1-barcodes.tsv")
peaks_primary <- read.table("CM618C1-S1-peaks.bed", sep = "\t")
counts_primary  <- Matrix::readMM("CM618C1-S1-matrix.mtx")
peaks_primary_gr <- GRanges(
  seqnames = peaks_primary$V1,
  ranges = IRanges(start = peaks_primary$V2, end = peaks_primary$V3)
)
peaknames_primary <- paste0(peaks_primary$V1, ":",peaks_primary$V2,"-",peaks_primary$V3)
colnames(counts_primary) <- barcodes_primary
rownames(counts_primary) <- peaknames_primary


fragpath_primary <- 'CM618C1-S1-fragments.tsv.gz'

rds_primary <- readRDS("CRC_CM618C2-S1Y2.rds")
metadata_primary <- rds_primary@meta.data

chrom_assay_primary <- CreateChromatinAssay(
  counts = counts_primary,
  sep = c(":", "-"),
  fragments = fragpath_primary,
  min.cells = 10,
  min.features = 200
)

primary_obj <- CreateSeuratObject(
  counts = chrom_assay_primary,
  assay = "peaks",
  meta.data = metadata_primary
)

primary_obj <- subset(
  primary_obj,
  cells = rds_primary$Original_barcode
)
cell_type_vector_names <- gsub(pattern = "CRC_CM618C2-S1Y2_",replacement = "",x = names(rds_primary$cell_type))
cell_type_vector <- rds_primary$cell_type
names(cell_type_vector) <- cell_type_vector_names
primary_obj$cell_type <- cell_type_vector
tumour_cells_primary <- names(which(primary_obj$cell_type=="Tumor"))
tumour_cells_met <- names(which(met_obj$cell_type=="Tumor"))

############ MERGE
met_obj <- RunTFIDF(met_obj)
met_obj <- FindTopFeatures(met_obj, min.cutoff = 'q0')
met_obj <- RunSVD(met_obj)
DepthCor(met_obj)

primary_obj <- RunTFIDF(primary_obj)
primary_obj <- FindTopFeatures(primary_obj, min.cutoff = 'q0')
primary_obj <- RunSVD(primary_obj)
DepthCor(primary_obj)


# first add dataset-identifying metadata
primary_obj$dataset <- "primary"
met_obj$dataset <- "metastasis"

# merge
crc.combined <- merge(met_obj, primary_obj)
tumour_cells_combined <- names(which(crc.combined$cell_type=="Tumor"))

# process the combined dataset
crc.combined <- FindTopFeatures(crc.combined, min.cutoff = 10)
crc.combined <- RunTFIDF(crc.combined)
crc.combined <- RunSVD(crc.combined)
crc.combined <- RunUMAP(crc.combined, reduction = "lsi", dims = 2:30)







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
Annotation(crc.combined) <- annotations



p1 <- DimPlot(crc.combined,split.by="dataset",group.by = "cell_type")
p4 <- DimPlot(crc.combined,split.by="dataset")
integrated <- RunHarmony(
  object = crc.combined,
  group.by.vars = 'dataset',
  reduction.use = 'lsi',
  # assay.use = 'peaks',
  project.dim = FALSE
)

integrated <- RunUMAP(integrated, dims = 2:30, reduction = 'harmony', reduction.name = 'harmony.umap')
integrated <- FindNeighbors(integrated, reduction = 'harmony', dims = 2:30)
integrated <- FindClusters(integrated, verbose = FALSE, algorithm = 3, resolution = 0.8)
p2 <-DimPlot(integrated, reduction = "harmony.umap", split.by="dataset",group.by = "cell_type")

markers <- FindAllMarkers(integrated,
                          assay="peaks",
                          min.pct = 0.2,
                          test.use = 'LR',
                          latent.vars = 'nCount_peaks',
                          only.pos = T, 
                          logfc.threshold = 0.25)


CoveragePlot(object = crc.combined,region = "ATRX",
             extend.upstream = 1000,extend.downstream = 1000,
             group.by = "cell_type",peaks = T,annotation = T)
