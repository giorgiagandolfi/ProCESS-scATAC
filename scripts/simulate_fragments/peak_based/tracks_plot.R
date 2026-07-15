library(ArchR)
library(ggplot2)
library(GenomicRanges)
library(dplyr)
addArchRThreads(threads = 16) 
addArchRGenome("hg38")



cell_epigenetic_df = read.table("/data/rds/DMP/UCEC/GENEVOD/ccolson/MECCA/AllCells/HallmarkPathwayAnalysis/HIPEC_P05/C4_Normalised/magic_gsva_bootstrap_zscoreParam_consensus_clusters.csv",header = T,sep=",")
archr_prj <- loadArchRProject(path = "/data/rds/DMP/UCEC/GENEVOD/ccolson/MECCA/AllCells/ArchRProjects/HIPEC_P05_AllCells/", force = FALSE, showLogo = TRUE)
peaks = readRDS("/data/rds/DMP/UCEC/GENEVOD/ccolson/MECCA/AllCells/ArchRProjects/HIPEC_P05_AllCells/PeakCalls/Type/Tumour-reproduciblePeaks.gr.rds")%>% as.data.frame()

myc_peaks <- peaks %>% 
  filter(nearestGene%in%c("CBS"))
  # filter(score>20)

myc_peaks_gr <- GRanges(
  seqnames = myc_peaks$seqnames,
  ranges = IRanges(start = myc_peaks$start, end = myc_peaks$end),
  strand = myc_peaks$strand
)


myc_peaks_not_filtered <- peaks %>% 
  filter(nearestGene%in%c("MYC","TGFA")) %>% 
  filter(score<=20)

myc_peaks_not_filtered_gr <- GRanges(
  seqnames = myc_peaks_not_filtered$seqnames,
  ranges = IRanges(start = myc_peaks_not_filtered$start, end = myc_peaks_not_filtered$end),
  strand = myc_peaks_not_filtered$strand
)

myc_grl <- GRangesList(
  filtered = myc_peaks_gr,
  all = myc_peaks_not_filtered_gr
)
pdf("myc_gene_p05.pdf")
archr_prj@cellColData$cell_id = rownames(archr_prj@cellColData)
cell_epigenetic_df = cell_epigenetic_df %>% mutate(
  cell_id = factor(cell_id, levels = rownames(archr_prj@cellColData))
) %>% 
  filter(!is.na(consensus_cluster)) %>% 
  tibble::column_to_rownames("cell_id")

consensus_cluster <- as.character(cell_epigenetic_df$consensus_cluster)
names(consensus_cluster) <- rownames(cell_epigenetic_df)
archr_prj <- addCellColData(
  ArchRProj = archr_prj,
  data = consensus_cluster,
  name = "consensus_cluster",
  cells = names(consensus_cluster),
  force = TRUE
)

pdf("myc_track.pdf")
p=plotBrowserTrack(
  ArchRProj = archr_prj, 
  groupBy = "consensus_cluster",
  geneSymbol = c("CBS"), 
  # useMatrix = "GeneScoreMatrix",
  upstream = 40000,
  downstream = 40000,
  features = myc_peaks_gr
)

grid.arrange(
  p$CBS,
  # p$TGFA,
  ncol = 1
)

dev.off()
