library(ArchR)
library(ComplexHeatmap)
set.seed(1)
addArchRThreads(threads = 16) 
addArchRGenome("hg38")

ArrowFiles <- createArrowFiles(
  inputFiles = 'analysis_bg/results_v7/fragment_calling/P05/P05.sorted.bed.gz',
  # inputFiles = "test_nextflow_pipeline/fragment_calling/sample_A/sample_A.chr.sorted.bed.gz",
  sampleNames = "P05",
  minFrags = 1,
  filterFrags = 1000,
  minTSS = 0,
  filterTSS = 4,
  addTileMat = T,
  addGeneScoreMat =T,excludeChr = c("chrY","chrX"),force = T
)

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "P05",
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)
cell_info = readRDS("0_process_simulations/sampled_cells_info_atac_epigenome_1.3.5_pat05.rds")
cell_info <- cell_info %>% 
  dplyr::filter(!is.na(sample)) %>% 
  # dplyr::mutate(class=paste(mutant,epistate,sep = "")) %>% 
  dplyr::mutate(sample='P05') %>% 
  dplyr::mutate(cell_id_correct=paste0(sample,"#",cell_id)) %>% 
  dplyr::filter(cell_id_correct%in%rownames(proj@cellColData))
cell_info <- cell_info[match(rownames(proj@cellColData), cell_info$cell_id_correct), ]


proj <- addCellColData(
  ArchRProj = proj,
  data = cell_info$epistate,
  cells = cell_info$cell_id_correct,
  name = "epigenetic_clade",
  force = TRUE
)

# doubScores <- addDoubletScores(
#   input = ArrowFiles,
#   k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
#   knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
#   LSIMethod = 1
# )

# proj <- addGeneScoreMatrix(proj,extendUpstream = 10,extendDownstream = 10)

proj_test <- getTestProject()

# Plot TSS
pt <- plotTSSEnrichment(proj)





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
plotFragmentSizes(ArchRProj = proj,groupBy = "Sample")

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
proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI")

p <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "epigenetic_clade", 
  geneSymbol = "MYC", 
  upstream = 10000,
  downstream = 10000
)
grid::grid.newpage()
grid::grid.draw(p$MYC)



proj@peakAnnotation


geneScores <- getMatrixFromProject(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix"
)



genes_hallmark_list = readRDS("/data/rds/DMP/UCEC/GENEVOD/ccolson/MECCA/AllCells/HallmarkPathwayAnalysis/HallmarkPathways.rds")
genes_hallmark <- data.frame(
  pathway = rep(names(genes_hallmark_list), lengths(genes_hallmark_list)),
  gene = unlist(genes_hallmark_list, use.names = FALSE)
)
mat <- assay(geneScores)
rownames(mat) <- rowData(geneScores)$name


selected_genes=genes_hallmark %>% 
  dplyr::filter(pathway%in%c('HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION','HALLMARK_MYC_TARGETS_V1')) %>% 
  dplyr::pull(gene) %>% 
  unique()
sampled_gene_prj=which(rownames(mat)%in%selected_genes)
mat_subsampled=mat[sampled_gene_prj,]


epigenetic_cluster_cols = c(
  "E1" = "olivedrab",
  "E2" = "forestgreen",
  "E3" = "darkgreen",
  "E4" ="palegreen1"
)


col_ann=HeatmapAnnotation(epigenetic_clade = proj@cellColData$`epigenetic clade`,col = list("epigenetic_clade" = epigenetic_cluster_cols))

ht=Heatmap(as.matrix(mat_subsampled),top_annotation = col_ann,show_column_names = F,
           show_column_dend = F,show_row_names = F,show_row_dend = F,
        col = colorRampPalette(c("white", "red"))(10))
pdf("simulted_gene_scores_archr.pdf")
draw(ht)
dev.off()
proj_chole = loadArchRProject('/data/rds/DMP/UCEC/GENEVOD/ccolson/MECCA/AllCells/ArchRProjects/HIPEC_P05_AllCells/')
geneScores_GT <- getMatrixFromProject(
  ArchRProj = proj_chole,
  useMatrix = "GeneScoreMatrix"
)

peakMatr_GT <- getMatrixFromProject(
  ArchRProj = proj_chole,
  useMatrix = "PeakMatrix"
)
peak_mat <- assay(peakMatr_GT)
rownames(peak_mat) <- rowData(peakMatr_GT)$name

peaks = readRDS("/data/rds/DMP/UCEC/GENEVOD/ccolson/MECCA/AllCells/ArchRProjects/HIPEC_P05_AllCells/PeakCalls/Type/Tumour-reproduciblePeaks.gr.rds")%>% as.data.frame()
peaks = peaks %>% 
  dplyr::mutate(peak=paste(seqnames,start,end,sep = "-"))
selected_peaks = readRDS("input_data_P05/data/final_version/peak_pathway_list.rds") %>% 
  dplyr::bind_rows() %>% dplyr::select(!pathway) %>% dplyr::distinct() %>% 
  dplyr::mutate(chr=paste0("chr",chr))

ggvenn::ggvenn(list("all_peaks"=peaks$peak,"selected_peaks"=selected_peaks$peak),auto_scale = T)


fragments =getFragmentsFromProject(ArchRProj = proj_chole)
available_samples <- names(fragments)
selected_peaks_gr = GRanges(
  seqnames = selected_peaks$chr,
  ranges = IRanges(start = selected_peaks$from, end = selected_peaks$to),
)
mapped_fragments_peaks <- list()
for (s in available_samples){
  fragments_sample <- fragments[[s]]
  overlaps_in_peak = findOverlaps(query = selected_peaks_gr,subject = fragments_sample,ignore.strand = TRUE)
  fragments_in_peaks <- fragments_sample[subjectHits(overlaps_in_peak)]
  fragments_sample=as.data.frame(fragments_sample)
  
  overlpaed_frags =fragments_sample[subjectHits(overlaps_in_peak),]
  overlpaed_peaks =selected_peaks[queryHits(overlaps_in_peak),]
  colnames(overlpaed_peaks) = paste0("peak_",colnames(overlpaed_peaks))
  colnames(overlpaed_frags) = paste0("fragment_",colnames(overlpaed_frags))
  joined = cbind(overlpaed_peaks,overlpaed_frags)
  joined$overlap_bp <- pmin(joined$fragment_end, joined$peak_to) -
    pmax(joined$fragment_start, joined$peak_from) + 1
  mapped_fragments_peaks[[s]] <- joined %>% 
    dplyr::mutate(out_peak_bp=abs(fragment_width-overlap_bp)) %>% 
    dplyr::mutate(pct_overlap=overlap_bp/fragment_width) %>% 
    dplyr::mutate(sample=s)
  
}
joined_all <- do.call("rbind",mapped_fragments_peaks)




joned_summary_statts = joined_all %>% 
  dplyr::group_by(peak_peak,sample) %>% 
  dplyr::mutate(n_fragments=dplyr::n()) %>% 
  dplyr::mutate(mean_overlap=mean(pct_overlap)) %>% 
  dplyr::mutate(n_cells= dplyr::n_distinct(fragment_RG)) %>% 
  dplyr::select(peak_peak,n_fragments,mean_overlap,n_cells,sample) %>% 
  dplyr::distinct() %>% 
  dplyr::mutate(n_fragments_per_cell=n_fragments/n_cells)


joined_all %>% dplyr::group_by(fragment_RG,sample) %>% 
  dplyr::summarise(tot_fragments_per_cell=dplyr::n()) %>% 
  ggplot(aes(x=tot_fragments_per_cell,fill=sample))+
  geom_boxplot(outliers = F)+theme_minimal()


joined_all %>% dplyr::group_by(fragment_RG,sample) %>% 
  dplyr::summarise(mean_overlap_per_cell=mean(pct_overlap)) %>% 
  ggplot(aes(x=mean_overlap_per_cell,fill=sample))+
  geom_boxplot(outliers = F)+theme_minimal()

joined_all %>% 
  dplyr::group_by(peak_peak,fragment_RG) %>% 
  dplyr::summarise(n_fragments_per_cell=dplyr::n()) %>% 
  ggplot(aes(x=n_fragments_per_cell))+
  geom_boxplot(outliers=F)




mat_gt <- assay(geneScores_GT)
rownames(mat_gt) <- rowData(geneScores_GT)$name

peakMat_GT <- getMatrixFromProject(
  ArchRProj = proj_chole,
  useMatrix = "PeakMatrix"
)
mat_peak_gt <- assay(peakMat_GT)

zeros_per_row <- (nrow(mat_peak_gt) - colSums(mat_peak_gt != 0))/nrow(mat_peak_gt)
mean(zeros_per_row)
library(fitdistrplus)
fit_lnorm <- fitdist(zeros_per_row, "lnorm")
fit_gamma <- fitdist(zeros_per_row, "gamma")
fit_weibull <- fitdist(zeros_per_row, "weibull")

summary(fit_lnorm)
summary(fit_gamma)
summary(fit_weibull)



gofstat(list(
  lognormal = fit_lnorm,
  gamma = fit_gamma,
  weibull = fit_weibull
))

plot(fit_lnorm)
plot(fit_weibull)
plot(fit_gamma)

fit_weibull$estimate


sampled_gene_prj=which(rownames(mat_gt)%in%selected_genes)
mat_gt_subsampled=mat_gt[sampled_gene_prj,]


Heatmap(as.matrix(mat_gt_subsampled),show_column_names = F,show_column_dend = F,show_row_names = F,show_row_dend = T,
        col = colorRampPalette(c("white", "red"))(10))


gene_score_gt <- as.data.frame.table(as.matrix(mat_gt_subsampled))
colnames(gene_score_gt) <- c("gene", "cell_id", "gene_score")

gene_score <- as.data.frame.table(as.matrix(mat_subsampled))
colnames(gene_score) <- c("gene", "cell_id", "gene_score")


gene_score_gt <- as.data.frame.table(as.matrix(mat_gt))
colnames(gene_score_gt) <- c("gene", "cell_id", "gene_score")

gene_score <- as.data.frame.table(as.matrix(mat))
colnames(gene_score) <- c("gene", "cell_id", "gene_score")

gene_score_gt <- gene_score_gt %>% 
  dplyr::group_by(gene) %>% 
  dplyr::summarise(mean_gene_score=mean(gene_score))


gene_score <- gene_score %>% 
  dplyr::group_by(gene) %>% 
  dplyr::summarise(mean_gene_score=mean(gene_score))

joined_gene_scores = dplyr::inner_join(gene_score,y = gene_score_gt,by='gene',suffix=c('_SIM','_GT'))
joined_gene_scores %>% 
  ggplot(aes(x=mean_gene_score_SIM,y=mean_gene_score_GT))+
  geom_point(alpha=0.2)+
  # geom_smooth()+
  geom_rug()+
  theme_minimal()+
  ggtitle(label = 'Correlation between ground truth and inferred gene scores')

cor(joined_gene_scores$mean_gene_score_SIM,y = joined_gene_scores$mean_gene_score_GT)



gene_score_gt_missing_values =
  gene_score_gt %>% dplyr::group_by(gene) %>% dplyr::summarise(
  n_zero = sum(gene_score == 0),
  n_nonzero = sum(gene_score != 0)) %>% 
  dplyr::mutate(pct_zeros=n_zero/(n_zero+n_nonzero))

gene_score_missing_values =
  gene_score %>% dplyr::group_by(gene) %>% dplyr::summarise(
    n_zero = sum(gene_score == 0),
    n_nonzero = sum(gene_score != 0)) %>% 
  dplyr::mutate(pct_zeros=n_zero/(n_zero+n_nonzero))


joined_gene_scores_missing_vals = dplyr::inner_join(gene_score_missing_values,y = gene_score_gt_missing_values,by='gene',suffix=c('_SIM','_GT'))
joined_gene_scores_missing_vals %>% 
  ggplot(aes(x=pct_zeros_SIM))+
  geom_boxplot()+
  theme_minimal()

joined_gene_scores_missing_vals %>% 
  ggplot(aes(x=pct_zeros_GT))+
  geom_boxplot()+
  theme_minimal()
library(tidyr)


joined_gene_scores_missing_vals %>% 
  dplyr::select(gene,pct_zeros_GT,pct_zeros_SIM) %>% 
  pivot_longer(
    cols = c(pct_zeros_GT, pct_zeros_SIM),
    names_to = "type",
    values_to = "pct_zeros"
  ) %>% 
  dplyr::mutate(type = dplyr::case_match(
    type,
    "pct_zeros_GT" ~ "Real",
    "pct_zeros_SIM" ~ "Simulated"
  )) %>% 
  ggplot(aes(x=pct_zeros,fill=type))+geom_boxplot()+
  theme_minimal()
