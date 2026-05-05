
library(Signac)
library(Seurat)
library(dplyr)
library(ggplot2)
#counts <- Read10X_h5(filename = "10k_pbmc_ATACv2_nextgem_Chromium_Controller_filtered_peak_bc_matrix.h5")
fragpath <- '10k_pbmc_ATACv2_nextgem_Chromium_Controller_fragments.tsv.gz'
#t_fr <- read.table(file = "CM663C1-T1-atac_fragments.tsv.gz",header = T,skip = 51)


# Define cells
# If you already have a list of cell barcodes to use you can skip this step
total_counts <- CountFragments(fragpath)
cutoff <- 1000 # Change this number depending on your dataset!
barcodes <- total_counts[total_counts$frequency_count > cutoff, ]$CB

# Create a fragment object
frags <- CreateFragmentObject(path = fragpath, cells = barcodes)

# First call peaks on the dataset
# If you already have a set of peaks you can skip this step
peaks <- CallPeaks(frags,
                   macs2.path = "/orfeo/cephfs/home/cdslab/ggandolfi/.cache/R/basilisk/1.18.0/zellkonverter/1.16.0/zellkonverterAnnDataEnv-0.10.9/bin/macs3")

# Quantify fragments in each peak
counts <- FeatureMatrix(fragments = frags, features = peaks, cells = barcodes)


getGEOSuppFiles(GEO = "GSE216175",fetch_files = T)

pbmc[['peaks']]
peaks.keep <- seqnames(granges(pbmc)) %in% GenomeInfoDb::standardChromosomes(granges(pbmc))
pbmc <- pbmc[as.vector(peaks.keep), ]
system("which python3")
system("which macs3")
system("python3 -m pip install --user MACS3")
