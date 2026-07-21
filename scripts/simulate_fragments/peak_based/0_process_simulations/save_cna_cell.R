library(ProCESS)
library(dplyr)

phylo_forest <- load_phylogenetic_forest("phylo_forest_atac_epigenome_1.3.5_pat05.sff")
single_cell_cnas = phylo_forest$get_cell_allelic_fragmentation() %>%
  mutate(karyotype=paste(major,minor,sep=":"))

saveRDS(object = single_cell_cnas,file = "cell_fragmentation.rds")