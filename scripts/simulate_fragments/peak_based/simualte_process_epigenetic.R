rm(list=ls())
library(ProCESS)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
# source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC/plot_muller.R")
# source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC/utils/DLP.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/scripts/utils/utils_plot.R")
setwd("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-scATAC/scripts/simulate_fragments/peak_based/")
sim_id <- "ATAC4"
set.seed(0)

mutant_cols <- c(
  "A+" = "goldenrod",
  "A-" = "magenta4",
  "B+" = "royalblue3",
  "B-" ="forestgreen"
)

sample_cols <- c(
  S1 = "lightskyblue1",
  S2 = "royalblue3"
)

sim <- TissueSimulation(name = sim_id,width = 500,height = 500,save_snapshots = F)

# set the delta time between two time series samples
sim$history_delta <- 10

# avoid drift
sim$death_activation_level <- 50

# sim$add_mutant(name = "A", growth_rates = 0.12, death_rates = 0.05)
sim$add_mutant(name = "A",
               epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
               growth_rates = c("+" = 0.2, "-" = 0.08),
               death_rates = c("+" = 0.05, "-" = 0.01))
sim$place_cell("A+", 250, 250)
sim$run_up_to_size("A+", 2000)
sim$run_up_to_size("A-", 2000)
sim$add_mutant(name = "B",
               epigenetic_rates = c("+-" = 0.7, "-+" = 0.04),
               growth_rates = c("+" = 0.4, "-" = 0.08),
               death_rates = c("+" = 0.05, "-" = 0.01))
sim$mutate_progeny(sim$choose_cell_in("A"), "B")

sim$run_up_to_size("B-", 2000)
plot_tissue(sim,color_map = mutant_cols)+
  facet_wrap(~species)

n_w <- n_h <- 10
ncells <- 0.5 * n_w * n_h

# Sampling ncells with random box sampling of boxes of size n_w x n_h
bbox1 <- sim$search_sample(c("A" = ncells), n_w, n_h)
bbox2 <- sim$search_sample(c("B" = ncells), n_w, n_h)

# plot the bounding box
s1 <- plot_tissue(sim, color_map = mutant_cols) +
  geom_rect(xmin = bbox1$lower_corner[1], xmax = bbox1$upper_corner[1],
            ymin = bbox1$lower_corner[2], ymax = bbox1$upper_corner[2],
            fill = NA, color = sample_cols["S1"])+
  geom_rect(xmin = bbox2$lower_corner[1], xmax = bbox2$upper_corner[1],
            ymin = bbox2$lower_corner[2], ymax = bbox2$upper_corner[2],
            fill = NA, color = sample_cols["S2"])
sim$sample_cells("S1", bbox1$lower_corner, bbox1$upper_corner)
sim$sample_cells("S2", bbox2$lower_corner, bbox2$upper_corner)
# sim$add_mutant(name = "B", growth_rates = 0.3, death_rates = 0.01)
# sim$mutate_progeny(sim$choose_cell_in("A"), "B")
# 
# sim$run_up_to_size("B", 2000)
# bbox <- sim$search_sample(c("B" = ncells*0.2,"A" = ncells*0.8), n_w, n_h)
# 
# 
# s1 <-plot_tissue(sim, color_map = mutant_cols) +
#   geom_rect(xmin = bbox$lower_corner[1], xmax = bbox$upper_corner[1],
#             ymin = bbox$lower_corner[2], ymax = bbox$upper_corner[2],
#             fill = NA, color = sample_cols["S1"])
# sim$sample_cells("S1", bbox$lower_corner, bbox$upper_corner)  
# 
# sim$add_mutant(name = "C", growth_rates = 0.4, death_rates = 0.01)
# sim$mutate_progeny(sim$choose_cell_in("A"), "C")
# sim$update_rates(name="B", c(growth = 0.150))
# 
# sim$run_up_to_size("C", 5000)
# bbox <- sim$search_sample(c("C" = ncells), n_w, n_h)
# 
# s2 <-plot_tissue(sim, color_map = mutant_cols) +
#   geom_rect(xmin = bbox$lower_corner[1], xmax = bbox$upper_corner[1],
#             ymin = bbox$lower_corner[2], ymax = bbox$upper_corner[2],
#             fill = NA, color = sample_cols["S2"])
# sim$sample_cells("S2", bbox$lower_corner, bbox$upper_corner)  
sample_forest <- sim$get_sample_forest()
plot_forest(forest = sample_forest,color_map = mutant_cols)

sample_forest$save("sample_forest_atac_epigenome.sff")
setwd("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC/")
m_engine <- MutationEngine(setup_code = "GRCh38",tumour_type = "COADREAD", context_sampling = 20)


mu_SNV = 1e-9
mu_CNA = 5e-12
mu_INDELs = 1e-9

CNA_Clone2 = ProCESS::CNA(type = "D", "5",
                          from = 107707518, len = 2e7)
CNA_Clone1 = ProCESS::CNA(type = "A", "12",
                          from = 10000000, len = 4e7)

## Drivers for the tumors
m_engine$add_mutant(mutant_name = "A",
                    passenger_rates = list("+" = c(SNV = mu_SNV, indel = mu_INDELs, CNA=mu_CNA),
                                           "-" = c(SNV = mu_SNV, indel = mu_INDELs, CNA=mu_CNA)),
                    drivers = list("KRAS G12D"))
m_engine$add_mutant(mutant_name = "B",
                    passenger_rates = list("+" = c(SNV = mu_SNV, indel = mu_INDELs, CNA=1e-11),
                                           "-" = c(SNV = mu_SNV, indel = mu_INDELs, CNA=mu_CNA)),
                    drivers = list("APC R1450*"))

# m_engine$add_mutant(mutant_name = "C",passenger_rates = c(SNV = mu_SNV, CNA = mu_CNA,indel=mu_INDELs),
#                     drivers = list("APC R1450*"))

m_engine$add_exposure(time = 0,coefficients = c(SBS1 = 0.15,SBS5 = 0.40,
                                                SBS18 = 0.15,SBS17b = 0.20,ID1 = 0.40,ID2 = 0.40,ID18=0.2,SBS88 = 0.10))
phylo_forest <- m_engine$place_mutations(sample_forest, num_of_preneoplatic_SNVs=800, num_of_preneoplatic_indels=200)
phylo_forest$save("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-scATAC/scripts/simulate_fragments/peak_based/phylo_forest_atac_epigenetic.sff")


# 
# 
# ggsave(filename = "muller_plot2.pdf",plot = muller_plot,width = 10,height = 4,dpi = 300)
# ggsave(filename = "phylo_forest_plot2.png",plot = plot_forest+theme_void(),width = 5,height = 5,dpi = 300)
# ggsave(filename = "s1_2.pdf",plot = s2,width = 4,height = 4,dpi = 300)
# ggsave(filename = "s2_2.pdf",plot = s3,width = 4,height = 4,dpi = 300)
# pdf(file = "heatmap.pdf",width = 5,height = 5)
# draw(plot_heatmap_chromtin)
# dev.off()
# 
# phylo_forest <- load_phylogenetic_forest("phylo_forest_atac1.sff")
# seq_results <- simulate_seq(phylo_forest, coverage = 10,chromosomes = c("5","12","17"),with_normal_sample = F)
