library(ProCESS)
library(dplyr)
library(ggplot2)
library(ggrepel)
setwd("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC/scripts/plot_muller.R")

########## TISSUE PLOTS
sim <-ProCESS::recover_simulation("ATAC2")
sample_forest <- ProCESS::load_sample_forest("sample_forest_atac2.sff")
phylo_forest <- ProCESS::load_phylogenetic_forest("phylo_forest_atac2.sff")
mutant_cols <- c(
  "Clone 1" = "goldenrod",
  "Clone 2" = "magenta4",
  "Clone 3" ="forestgreen"
)
plot_muller <- plot_Muller(sim,color_map = mutant_cols)
s1 <- plot_tissue(simulation = sim,at_sample = "S1",color_map = mutant_cols)
s2 <- plot_tissue(simulation = sim,at_sample = "S2",color_map = mutant_cols)
time_series <- plot_timeseries(sim,color_map = mutant_cols)

sample_composition <- samples_table(sim=sim,use_snapshot = F, sample_forest=sample_forest)
sample_composition_long <- sample_composition %>%
  select(Sample_ID, "Clone 2", "Clone 3") %>%
  pivot_longer(cols = c("Clone 2", "Clone 3"),
               names_to = "mutation",
               values_to = "count")

sample_comp_pltos <- ggplot(sample_composition_long,
       aes(x = "", y = count, fill = mutation)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  facet_wrap(~Sample_ID,scales="free") +
  theme_void() +
  labs(fill = "Mutant")+
  scale_fill_manual(values=mutant_cols)

plot_forest <- plot_forest(sample_forest,color_map = mutant_cols) %>%
  my_annotate_forest(phylo_forest,add_driver_label =F,drivers = F,color_map = mutant_cols)

#ggsave(filename = "muller_plot2.pdf",plot = plot_muller,width = 6,height = 3,dpi = 300)
dir.create(path =  "/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC/mets_plot",recursive = T)
ggsave(filename = "/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC/mets_plot/phylo_forest_plot2.png",plot = plot_forest,width = 5,height = 8,dpi = 300)
ggsave(filename = "/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC/mets_plot/s1_2.png",plot = s1,width = 4,height = 4,dpi = 300)
ggsave(filename = "/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC/mets_plot/s2_2.png",plot = s2,width = 4,height = 4,dpi = 300)
ggsave(filename = "/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC/mets_plot/time_series.pdf",plot = time_series,width = 4,height = 4,dpi = 300)
ggsave(filename = "/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC/mets_plot/sample_comp_plots.pdf",plot = sample_comp_pltos,width = 4,height = 4,dpi = 300)

seq_results <- readRDS("seq_results3.rds")
seq_results_muts <- lapply(seq_along(1:length(seq_results)), function(i){
  muts <- seq_results[[i]]$mutations
}) %>% bind_rows()

sample_names <- strsplit(colnames(seq_results_muts)[grepl(".VAF",
                                                          colnames(seq_results_muts),
                                                          fixed = TRUE)],
                         ".VAF") %>% unlist()

# Process each sample separately to create a list of data frames
seq_df <- lapply(sample_names, function(sn) {
  # Select relevant columns for the current sample
  cc <- c("chr", "from", "ref", "alt", "causes", "classes",
          colnames(seq_results_muts)[grepl(paste0(sn, "."),
                                           colnames(seq_results_muts), fixed = TRUE)])
  
  # Rename columns and add sample_name column
  seq_results_muts[, cc] %>%
    `colnames<-`(c("chr", "from", "ref", "alt", "causes", "classes",
                   "occurrences", "coverage", "VAF")) %>%
    dplyr::mutate(sample_name = sn)
}) %>% do.call("bind_rows", .)

# Rename and reorder columns
seq_results_muts_long <- seq_df %>%
  dplyr::rename(DP = coverage,
                NV = occurrences) %>%
  dplyr::mutate(to = from)



phylo_forest <- load_phylogenetic_forest("phylo_forest_atac3.sff")

cna_seg <- lapply(sample_names, FUN = function(s){
  cna_data <- phylo_forest$get_bulk_allelic_fragmentation(s)
  cna_data %>%
    dplyr::mutate(CN_type = ifelse(ratio < 0.9 & ratio > 0.1, 'sub-clonal', 'clonal'),
                  CN = paste(major, minor, sep = ':'),
                  seg_id = paste(chr,begin,end, sep = ':'),
                  sample = s)
})
names(cna_seg) <- sample_names
print('Done')

print('Map mutations to cnas')
muts_cn <- list()
for (sample in sample_names){
  print(sample)
  tmp_cna <- cna_seg[[sample]]
  
  muts_cn[[sample]] <- seq_results_muts_long %>%
    filter(sample_name==sample) %>% 
    filter(classes != "germinal") %>% 
    inner_join(tmp_cna, by = c("chr"), relationship = "many-to-many") %>%
    filter(from >= begin & to <= end)
}
muts_cn_df <- do.call("rbind",muts_cn)

muts_cn_df %>% 
  filter(VAF>=0.1) %>% 
  filter(CN=="1:1") %>% 
  ggplot(aes(x=VAF))+
  geom_histogram(binwidth = 0.01)+
  facet_wrap(~sample)

drivers_muts <- phylo_forest$get_driver_mutations() %>% 
  rename(from=start)
seq_res_somatic <- seq_results_muts %>% 
  filter(classes!="germinal") %>% 
  left_join(drivers_muts)
multivariate_vaf <- seq_res_somatic %>% 
  ggplot(aes(x=S1.VAF,y=S2.VAF))+
  geom_point(color="grey",alpha=0.2)+
  # gghighlight::gghighlight(classes == "driver")+
  geom_label_repel(
    data = ~ filter(.x, classes == "driver") %>% unique(),
    aes(label = code,fill=causes),
    # arrow = arrow(
    #   length = unit(2, "mm"), type = "closed", ends = "last"
    # ),
    size = 3,
    alpha=0.5,
    max.overlaps = Inf,
    min.segment.length = 0,
    box.padding = 0.9,
    point.padding = 0.5
  )+
  scale_fill_manual(values=mutant_cols)+
  theme_minimal()+
  labs(
    color = "Clone",
    fill = "Clone")+
  xlim(0,1)+
  ylim(0,1)

multivariate_vaf
marginal_plot_s1 <- seq_res_somatic %>% 
  filter(S1.VAF >= 0.05) %>% 
  ggplot(aes(x = S1.VAF)) +
  geom_histogram(binwidth = 0.01,fill="grey") +
  geom_vline(
    data = ~ filter(.x, classes == "driver"),
    aes(xintercept = S1.VAF,color=causes),
    linetype = "dashed",
  ) +
  # geom_label_repel(
  #   data = ~ filter(.x, classes == "driver"),
  #   aes(
  #     x = S1.VAF,
  #     y = Inf,
  #     label = code,
  #     fill=causes
  #   ),
  #   inherit.aes = FALSE,
  #   direction = "y",
  #   nudge_y = 5,
  #   vjust = 0,
  #   alpha=0.5,
  #   size = 3
  # ) +
  scale_color_manual(values=mutant_cols)+
  scale_fill_manual(values=mutant_cols)+
  theme_minimal()+
  labs(
    color = "Clone",
    fill = "Clone",
    y="Mutation count"
    # linetype = "Driver status"
  )

marginal_plot_s2 <- seq_res_somatic %>% 
  filter(S2.VAF >= 0.05) %>% 
  ggplot(aes(x = S2.VAF)) +
  geom_histogram(binwidth = 0.01,fill="grey") +
  geom_vline(
    data = ~ filter(.x, classes == "driver"),
    aes(xintercept = S2.VAF,color=causes),
    linetype = "dashed",
  ) +
  # geom_label_repel(
  #   data = ~ filter(.x, classes == "driver"),
  #   aes(
  #     x = S2.VAF,
  #     y = Inf,
  #     label = code,
  #     fill=causes
  #   ),
  #   inherit.aes = FALSE,
  #   direction = "y",
  #   nudge_y = 5,
  #   vjust = 0,
  #   alpha=0.5,
  #   size = 3
  # ) +
  scale_color_manual(values=mutant_cols)+
  scale_fill_manual(values=mutant_cols)+
  theme_minimal()+
  labs(
    color = "Clone",
    fill = "Clone",
    y="Mutation count"
    # linetype = "Driver status"
  )
wrap_plots(list(marginal_plot_s2,marginal_plot_s1,multivariate_vaf),
           design = "AACC\nBBCC",guides="collect") & theme(legend.position = "bottom")
ggsave(filename = "mets_plot/vaf_plots.pdf",width = 8,height = 4,dpi = 300)
