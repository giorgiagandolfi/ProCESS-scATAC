library(Seurat)
library(Signac)
library(Rsamtools)
library(GenomicRanges)
setwd("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC/")
## get data from the primary tumour
barcodes_primary <- readLines("CM618C1-S1-barcodes.tsv")
peaks_primary <- read.table("CM618C1-S1-peaks.bed", sep="\t")
peaks_primary_gr <- GRanges(
  seqnames = peaks_primary$V1,
  ranges = IRanges(start = peaks_primary$V2, end = peaks_primary$V3)
)
peaknames_primary <- paste(peaks_primary$V1, peaks_primary$V2, peaks_primary$V3, sep="-")
fragpath_primary <- 'CM618C1-S1-fragments.tsv.gz'
process_primary <- readRDS("CRC_CM618C2-S1Y2.rds")
barcodes_process_primary <- process_primary$Original_barcode

## get data from the met tumour
barcodes_met <- readLines("CM618C2-T1-barcodes.tsv")
peaks_met <- read.table("CM618C2-T1-peaks.bed", sep="\t")
peaks_met_gr <- GRanges(
  seqnames = peaks_met$V1,
  ranges = IRanges(start = peaks_met$V2, end = peaks_met$V3)
)
peaknames_met <- paste(peaks_met$V1, peaks_met$V2, peaks_met$V3, sep="-")
fragpath_met <- 'CM618C2-T1-fragments.tsv.gz'


##### combine the peaks
combined.peaks <- reduce(x = c(peaks_primary_gr,peaks_met_gr))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

# processed_data <- readRDS("CRC_CM268C1-T1.rds")
# cells <- sub(".*_", "", colnames(processed_data))
# Define cells
# If you already have a list of cell barcodes to use you can skip this step
# total_counts <- CountFragments(fragpath)
# cutoff <- 1000 # Change this number depending on your dataset!
# barcodes <- total_counts[total_counts$frequency_count > cutoff, ]$CB

# Create a fragment object
frags_mets <- CreateFragmentObject(path = fragpath_met, cells = barcodes_met)
frags_prim <- CreateFragmentObject(path = fragpath_primary, cells = barcodes_process_primary)

primary.counts <- FeatureMatrix(
  fragments = frags_prim,
  features = combined.peaks,
  cells = barcodes_process_primary
)

met.counts <- FeatureMatrix(
  fragments = frags_mets,
  features = combined.peaks,
  cells = barcodes_met
)



primary_assay <- CreateChromatinAssay(primary.counts, fragments = frags_prim)
primary_obj <- CreateSeuratObject(primary_assay, assay = "ATAC")

met_assay <- CreateChromatinAssay(met.counts, fragments = frags_mets)
met_obj <- CreateSeuratObject(met_assay, assay = "ATAC")



met_obj$dataset <- 'met'
primary_obj$dataset <- 'primary'


# merge all datasets, adding a cell ID to make sure cell names are unique
combined <- merge(
  x = primary_obj,
  y = met_obj,
  add.cell.ids = c("primary", "metastasis")
)
combined[["ATAC"]]


CoveragePlot(
  object = combined,
  group.by = 'dataset',peaks = T,
  region = "chr6-167855683-167856183",extend.upstream = 4000,extend.downstream = 4000
)



# First call peaks on the dataset
# If you already have a set of peaks you can skip this step
peaks <- CallPeaks(frags,macs2.path = "/orfeo/cephfs/home/cdslab/ggandolfi/.cache/R/basilisk/1.18.0/zellkonverter/1.16.0/zellkonverterAnnDataEnv-0.10.9/bin/macs3")

# Quantify fragments in each peak
counts <- FeatureMatrix(fragments = frags, features = peaks_gr, cells = barcodes)



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
)

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
Annotation(pbmc) <- annotations





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


regions_highlight <- subsetByOverlaps(StringToGRanges(open_cd4naive), LookupGeneCoords(pbmc, "CD4"))

CoveragePlot(
  object = pbmc,
  region = "SNAI1",
  # region.highlight = regions_highlight,
  extend.upstream = 1000,
  extend.downstream = 1000
)

pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = pbmc, group.by = 'nucleosome_group')
