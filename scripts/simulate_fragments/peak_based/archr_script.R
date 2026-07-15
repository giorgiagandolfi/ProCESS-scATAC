library(ArchR)
set.seed(1)
addArchRThreads(threads = 16) 
addArchRGenome("hg38")

ArrowFiles <- createArrowFiles(
  inputFiles = 'analysis_bg/results_v5/fragment_calling/sample_A/sample_A.fragments.tsv.gz',
  # inputFiles = "test_nextflow_pipeline/fragment_calling/sample_A/sample_A.chr.sorted.bed.gz",
  sampleNames = "sample_A",
  minFrags = 1,
  filterFrags = 1,
  minTSS = 0,
  filterTSS = 0,
  addTileMat = T,
  addGeneScoreMat =T,force = F,excludeChr = c("chrY","chrX")
)

# doubScores <- addDoubletScores(
#   input = ArrowFiles,
#   k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
#   knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
#   LSIMethod = 1
# )
proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "sample_A",
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)
# proj <- addGeneScoreMatrix(proj,extendUpstream = 10,extendDownstream = 10)




p1 <- plotGroups(
  ArchRProj = proj, 
  name = "TSSEnrichment",
  plotAs = "ridges",
  baseSize = 10
)
p3 <- plotGroups(
  ArchRProj = proj, 
  name = "log10(nFrags)",
  plotAs = "ridges",
  baseSize = 10
)
p4 <- plotTSSEnrichment(ArchRProj = proj)


proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30
)

projHeme2 <- addClusters(
  input = proj,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8
)


p <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "Sample", 
  geneSymbol = "EPHB2", 
  upstream = 50000,
  downstream = 50000
)
grid::grid.newpage()
grid::grid.draw(p$EPHB2)



proj@peakAnnotation


geneScores <- getMatrixFromProject(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix"
)

mat <- assay(geneScores)
rownames(mat) <- rowData(geneScores)$name
