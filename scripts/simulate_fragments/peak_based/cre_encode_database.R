library(AnnotationDbi)
library(org.Hs.eg.db)
genes_hallmark_list = readRDS("/data/rds/DMP/UCEC/GENEVOD/ccolson/MECCA/AllCells/HallmarkPathwayAnalysis/HallmarkPathways.rds")
genes_hallmark <- data.frame(
  pathway = rep(names(genes_hallmark_list), lengths(genes_hallmark_list)),
  gene = unlist(genes_hallmark_list, use.names = FALSE)
)
encode_genes <- read.table(file = 'GRCh38-Closest-Genes-PC.tsv')
colnames(encode_genes) <- c('ccre_chrom','ccre_start','ccre_end','ccre_id1','ccre_id2','ccre_type',
                            'gene_chrom','gene_start','gene_end','gene_tranID','gene_score',
                            'gene_strand','gene_id')
encode_genes$ensembl_id <- sub(
  "\\..*",
  "",
  encode_genes$gene_id
)

# convert ENSG -> SYMBOL
encode_genes$gene_symbol <- mapIds(
  org.Hs.eg.db,
  keys = encode_genes$ensembl_id,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)
encode_genes<- encode_genes %>% 
  mutate(distance_from_gene_start=abs(ccre_end-gene_start)) %>% 
  mutate(distance_from_gene_end=abs(ccre_start-gene_end)) %>% 
  filter(ccre_type%in%c("pELS","dELS","PLS"))


encode_genes <- encode_genes %>%
  filter(
    ccre_type %in% c("PLS","pELS","dELS")
  ) %>%
  mutate(
    distance = case_when(
      ccre_end < gene_start ~ gene_start - ccre_end,
      ccre_start > gene_end ~ ccre_start - gene_end,
      TRUE ~ 0
    )
  ) %>% 
  filter(distance < 100000) %>% 
  mutate(ccre_width=ccre_end-ccre_start)


encode_genes <- encode_genes %>%
  mutate(
    type_weight = case_when(
      ccre_type=="PLS"  ~ 1,
      ccre_type=="pELS" ~ 0.7,
      ccre_type=="dELS" ~ 0.4
    ),
    distance_weight = exp(-distance/50000),
    ccre_weight = type_weight * distance_weight
  )

notch_genes_reg_region <- encode_genes %>% 
  filter(gene_symbol%in%genes_hallmark_list$HALLMARK_NOTCH_SIGNALING)

notch_genes_reg_region %>% 
  filter(gene_symbol=="WNT2") %>% 
  ggplot(aes(x=ccre_id2,y=ccre_weight,fill=ccre_type))+
  geom_col()


glibrary(GenomicRanges)

# genes
gene_gr <- GRanges(
  seqnames = "chr8",
  ranges = IRanges(
    start = 127735434,
    end   = 127742951
  ),
  gene="MYC"
)

# cCREs
ccre <- read.table(file = 'GRCh38-cCREs.bed',header = F)
colnames(ccre) <- c('chr','start','end','v1','v2','v3')

ccre_gr <- GRanges(
  seqnames = ccre$chr,
  ranges=IRanges(
    start   = ccre$start,
    end   = ccre$end
  ),type=ccre$v3,
  v1=ccre$v1,
  v2=ccre$v2
)

hits <- distanceToNearest(
  ccre_gr,
  gene_gr
)

result <- data.frame(
  ccre=ccre_gr[queryHits(hits)],
  gene=gene_gr[subjectHits(hits)],
  distance=mcols(hits)$distance
)
result

