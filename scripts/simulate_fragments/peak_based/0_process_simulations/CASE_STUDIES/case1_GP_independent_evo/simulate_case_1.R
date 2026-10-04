rm(list=ls())
library(ProCESS)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggalluvial)
library(ComplexHeatmap)
library(patchwork)
setwd("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project/ProCESS-scATAC/scripts/simulate_fragments/peak_based/0_process_simulations/CASE_STUDIES/case1_GP_independent_evo/")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project/ProCESS-scATAC/scripts/simulate_fragments/peak_based/utils.R")
source("utils_case1.R")

set.seed(0)
dir.create("results", showWarnings = FALSE)
dir.create(file.path("results", "figures"), showWarnings = FALSE, recursive = TRUE)



sim <- TissueSimulation(epigenetic_states = c("P1", "P2","P3"),
                        width = 1000, height = 1000,save_directory = F,name = "case1")

rates_G1 <- list(
  P1 = list(duplication = 0.8, death = 0.1, P2 = 0.2,P3 = 0.6),
  P2 = list(duplication = 0.8, death = 0.1, P1 = 0.3,P3=0.6),
  P3 = list(duplication = 0.8, death = 0.1, P1 = 0.3,P2 = 0.2)
)


# add a mutant "A" and set its species rates
rates_G2 <- list(
  P1 = list(duplication = 1, death = 0.1,P2 = 0.2,P3 = 0.6),
  P2 = list(duplication = 1, death = 0.1, P1 = 0.3,P3=0.6),
  P3 = list(duplication = 1, death = 0.1, P1 = 0.3,P2=0.2)
)


rates_G3 <- list(
  P1 = list(duplication = 1.5, death = 0.1,P2 = 0.2,P3 = 0.6),
  P2 = list(duplication = 1.5, death = 0.1, P1 = 0.3,P3=0.6),
  P3 = list(duplication = 1.5, death = 0.1, P1 = 0.3,P2=0.2)
)

q_G1=get_epigenetic_Q(epigenetic_rates = rates_G1,epistate_name = "G1")
q_G2=get_epigenetic_Q(epigenetic_rates = rates_G2,epistate_name = "G2")
q_G3=get_epigenetic_Q(epigenetic_rates = rates_G3,epistate_name = "G3")


sim$add_mutant("G1", rates_G1)
sim$place_cell("G1[P1]", 500, 500)

# let the simulation evolve until the species "A[E2]" has less than 10 cells
sim$run_up_to_time(20)


sim$add_mutant("G2",rates_G2)
starting_cell=sim$choose_border_cell_in("G1")

sim$mutate_progeny(starting_cell,"G2")
sim$set_rates(list("G1" = list(duplication = 0.6)))
sim$run_up_to_time(40)


# Sampling ncells with random box sampling of boxes of size n_w x n_h

sim$add_mutant("G3",rates_G3)
starting_cell=sim$choose_border_cell_in("G2")

sim$mutate_progeny(starting_cell,"G3")
sim$set_rates(list("G2" = list(duplication = 0.6, death = 0.15)))
sim$set_rates(list("G1" = list(duplication = 0.4, death = 0.2)))
sim$run_up_to_time(80)

########## sample
bbox_width=40
bbox1_p <- c(350,430)
bbox1_q <- bbox1_p + bbox_width
sim$sample_cells("S3", bbox1_p, bbox1_q)
sample_forest <- sim$get_sample_forest()


#### plot report
sample_summary <- summarise_sample(sample_forest)
colors <- make_case1_color_maps()
sankey_epi <- plot_epistate_alluvial(sample_summary$epi_by_mutant, colors$epi)
sankey_gen <- plot_mutant_alluvial(sample_summary$mutant_by_epi, colors$gen)

p_forest_epi <- plot_forest(sample_forest, color_map=colors$epistate) +
  theme(legend.position="none") + ggtitle("Sample forest: phenotype")
p_forest_gen <- plot_forest(sample_forest, color_map=colors$mutant) +
  theme(legend.position="none") + ggtitle("Sample forest: genetic clone")


p_evo <- wrap_plots(list=list(p_forest_gen,p_forest_epi,sankey_gen,sankey_epi),design = "aabb\naabb\nccdd") +
  plot_annotation(tag_levels = "a")


p_tran_mat <- wrap_plots(list=list(q_G1$plot_q,q_G2$plot_q,q_G3$plot_q))+
  plot_annotation(tag_levels = "a")

ggsave("results/figures/01_gen_phen_evolution.pdf", p_evo, width=12, height=8, dpi=300)
ggsave("results/figures/02_transition_mat.pdf", p_tran_mat, width=12, height=4, dpi=300)
sample_forest$save("results/sample_forest.sff")





############## MUTATION ENGINE #################
dir.create("process_references_v1.3.5")
setwd("process_references_v1.3.5")

m_engine <- MutationEngine(setup_code = "GRCh38",tumour_type = "COADREAD", context_sampling = 20,
                           germline_subject = "NA20514",
                           COSMIC_account = list("email"="giorgia.gandolfi@phd.units.it","password"="2*db!XQ4sgQ!dbg"))


mu_SNV = 1e-9
mu_CNA = 0
mu_INDELs = 1e-9

CNA_Clone0_1 = ProCESS::CNA(type = "D", "11",
                            from = 48500001, len = 1e7)
CNA_Clone0_2  = ProCESS::CNA(type = "D", "5",
                             from = 1000001, len = 109499999,allele = 0)
CNA_Clone0_3 = ProCESS::CNA(type = "D","17",
                            from = 500001, len = 42999999)

# chrX      5 107000001 150500000       0.25         1 chrX:107000001:150500000 43499999 D

CNA_CloneA = ProCESS::CNA(type = "D", "X",
                          from = 107000001, len = 43499999)



CNA_CloneB_1 = ProCESS::CNA(type = "A", "8",
                            from = 127118340, len = 1e7,src_allele = 1)


## Drivers for the tumors
m_engine$add_mutant(mutant_name = "G1",
                    passenger_rates = list("E1" = c(SNV = mu_SNV, indel = mu_INDELs,CNA=mu_CNA),
                                           "E2" = c(SNV = mu_SNV, indel = mu_INDELs,CNA=mu_CNA),
                                           "E3" = c(SNV = mu_SNV, indel = mu_INDELs,CNA=mu_CNA),
                                           "E4" = c(SNV = mu_SNV, indel = mu_INDELs,CNA=mu_CNA)),
                    drivers = list("APC Q1294Gfs*6","TP53 R175H",CNA_Clone0_1,CNA_Clone0_2,CNA_Clone0_3))
m_engine$add_mutant(mutant_name = "G2",
                    passenger_rates = list("E1" = c(SNV = mu_SNV, indel = mu_INDELs,CNA=mu_CNA),
                                           "E2" = c(SNV = mu_SNV, indel = mu_INDELs,CNA=mu_CNA),
                                           "E3" = c(SNV = mu_SNV, indel = mu_INDELs,CNA=mu_CNA),
                                           "E4" = c(SNV = mu_SNV, indel = mu_INDELs,CNA=mu_CNA)),
                    drivers = list("KRAS Q22K",CNA_CloneA))

m_engine$add_mutant(mutant_name = "G3",
                    passenger_rates = list("E1" = c(SNV = mu_SNV, indel = mu_INDELs,CNA=mu_CNA),
                                           "E2" = c(SNV = mu_SNV, indel = mu_INDELs,CNA=mu_CNA),
                                           "E3" = c(SNV = mu_SNV, indel = mu_INDELs,CNA=mu_CNA),
                                           "E4" = c(SNV = mu_SNV, indel = mu_INDELs,CNA=mu_CNA)),
                    drivers = list("GNAS R844H",CNA_CloneB_1))

m_engine$add_exposure(time = 0,coefficients = c(SBS1 = 0.15,SBS5 = 0.40,
                                                SBS18 = 0.15,SBS17b = 0.20,ID1 = 0.40,ID2 = 0.40,ID18=0.2,SBS88 = 0.10))
phylo_forest <- m_engine$place_mutations(sample_forest, num_of_preneoplatic_SNVs=800, num_of_preneoplatic_indels=200)
plot_forest(phylo_forest,color_map = color_map_gen) %>%
  annotate_forest(phylo_forest,drivers = T,add_driver_label = T)
phylo_forest$save("../phylo_forest_atac_case_1.sff")

