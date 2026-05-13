rm(list = ls())
gc()

library(Seurat)
library(Signac)
# library(qs)
# library(ggpubr)
library(dplyr)

# setwd("/shared/ifbstor1/projects/denovo_signatures/PJ_MECCA/Script")

# define directory-----------

# 2) load the data -------
# dir_RNA <- file.path(dir_Terekhanova, "0_Processed_data_by_author", "snRNAseq")
# dfs <- list.files(dir_RNA)[!grepl("snRNA_L4", list.files(dir_RNA))]
dir_ATAC <- file.path("/data/rds/DMP/UCEC/GENEVOD/ggandolfi/HTAN_data/snATACseq/")
dfs <- list.files(dir_ATAC)[!grepl("snRNA_L4", list.files(dir_ATAC))]
dfs <- c("CRC_CM618C1-T1Y2.rds","CRC_CM618C2-S1Y2.rds")

seurat_obj_list <- list()
for (dat in dfs) {
  print(dat)
  tmp <- readRDS(file.path(dir_ATAC, dat))
  seurat_obj_list[[dat]] = tmp
}

colnames(tmp)


# keep previous UMAP provided by the authors-------- 
seurat_obj_list_rm_frag <- seurat_obj_list
# DefaultAssay(seurat_obj_list[[1]])
# "pancan"

# remove the fragment index
for (i in seq_along(seurat_obj_list_rm_frag)) {
  Fragments(seurat_obj_list_rm_frag[[i]][["pancan"]]) <- NULL
}

seurat_obj_list_rm_frag <-list(process_met,process_primary)
seurat_obj <- merge(seurat_obj_list_rm_frag[[1]], y = seurat_obj_list_rm_frag[-1],
                    add.cell.ids = NULL, merge.dr=TRUE,
                    project = "Terekhanova")


seurat_obj@meta.data <- seurat_obj@meta.data %>% 
  dplyr::mutate(sample_type=case_when(Piece_ID=="CRC_CM618C1-T1Y2"~"met",
                                      TRUE~"primary"))

saveRDS(object = seurat_obj,file="merged_seurat_obj.rds")
tumour_obj <- subset(
  seurat_obj,
  subset = cell_type == "Tumor"
)
Idents(tumour_obj) <- "sample_type"
table(Idents(tumour_obj))
my_cols <- c("darksalmon","forestgreen")
names(my_cols)<-c("primary","met")
cells_tumor <- seurat_obj@meta.data %>% filter(cell_type=="Tumor") %>% rownames()

p_umap <- DimPlot(seurat_obj, reduction = "umap",  
                  group.by = c("cell_type","sample_type"),
                  label = F, label.size = 3) #,cols = my_cols)
my_cols_cell_types <- c(
  "Tumor"                     = "#D73027",
  "Endothelial"               = "#1F78B4",
  "Hepatocytes"               = "#33A02C",
  "Fibroblasts"               = "#FB9A99",
  "Macrophages"               = "#FF7F00",
  "T-cells"                   = "#6A3D9A",
  "Cholangiocytes"            = "#B2DF8A",
  "Low quality"               = "#BDBDBD",
  "B-cells"                   = "#A6CEE3",
  "Plasma"                    = "#CAB2D6",
  "Distal Enterocytes"        = "#FDBF6F",
  "Goblet"                    = "#8DD3C7",
  "Distal TA"                 = "#FFFFB3",
  "Distal Mature Enterocytes" = "#BC80BD",
  "Distal Absorptive"         = "#CCEBC5",
  "Distal Stem Cells"         = "#FFED6F",
  "Other_doublets"            = "#636363"
)
p_umap <- DimPlot(seurat_obj, reduction = "umap",  
                  group.by = c("cell_type",'sample_type'),
                  
                  label = F, label.size = 3,cols = c(my_cols_cell_types,my_cols))

dar <- FindMarkers(
  object = tumour_obj,
  ident.1 = "met",
  ident.2 = "primary",
  test.use = "LR",
  latent.vars = "nCount_pancan",
  min.pct = 0.05,
  logfc.threshold = 0.1
)
saveRDS(object = dar,file = "dar_df_tumourcells.rds")
FeaturePlot(
  object = seurat_obj,
  features = "chr8-102653727-102654227",
  reduction = "umap",
)
dar = readRDS("dar_df_tumourcells.rds")
dar_top <- dar %>% head(10) %>% rownames()

gene.activities <- GeneActivity(tumour_obj)


peaks <- c(
  "chr20-62320630-62321130",
  "chr20-62321282-62321782"
)

# split into components
parts <- do.call(rbind, strsplit(peaks, "-"))

gr <- GRanges(
  seqnames = parts[,1],
  ranges = IRanges(
    start = as.numeric(parts[,2]),
    end   = as.numeric(parts[,3])
  )
)
p=CoveragePlot(
  object = seurat_obj,
  annotation = TRUE,
  group.by = "cell_type",
  peaks = TRUE,
  region = c(dar_top[1],'chr20-62320630-62321782'),
  # region.highlight = gr,
  extend.upstream = 1000,
  extend.downstream = 1000,
)
p= p &scale_fill_manual(values=my_cols_cell_types)

roi_matrix <- FeatureMatrix(
  fragments = Fragments(tumour_obj),
  features = gr,
  cells = colnames(tumour_obj)
)
tumour_obj$roi_accessibility <- as.numeric(roi_matrix[1, ])

ggsave(p, file = "selected_peaks.pdf", width = 10, height = 3)

library(AnnotationHub)
ah <- AnnotationHub()

# Search for the Ensembl 98 EnsDb for Homo sapiens on AnnotationHub
query(ah, "EnsDb.Hsapiens.v98")
ensdb_v98 <- ah[["AH75011"]]
seurat_obj = readRDS("merged_seurat_obj.rds")

annotations <- GetGRangesFromEnsDb(ensdb = ensdb_v98)

# change to UCSC style since the data was mapped to hg38
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"
Annotation(seurat_obj) <- annotations

gene.activities <- GeneActivity(tumour_obj)
tumour_obj[['RNA']] <- CreateAssayObject(counts = gene.activities)
tumour_obj <- NormalizeData(
  object = tumour_obj,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(tumour_obj$nCount_RNA)
)

DefaultAssay(pbmc) <- 'RNA'

FeaturePlot(
  object = pbmc,
  features = c('MS4A1', 'CD3D', 'LEF1', 'NKG7', 'TREM1', 'LYZ'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)

