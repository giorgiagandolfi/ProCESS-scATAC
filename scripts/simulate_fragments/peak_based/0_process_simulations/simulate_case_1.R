rm(list=ls())
library(ProCESS)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggalluvial)
library(ComplexHeatmap)
setwd("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project/ProCESS-scATAC/scripts/simulate_fragments/peak_based/0_process_simulations")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project/ProCESS-scATAC/scripts/simulate_fragments/peak_based/utils.R")
set.seed(0)
sim <- TissueSimulation(epigenetic_states = c("E1", "E2","E3"),
                        width = 1000, height = 1000,save_snapshots = F,name = "case1")

rates_genetic_clone_0 <- list(
  E1 = list(duplication = 0.8, death = 0.1, E2 = 0.2,E3 = 0.6),
  E2 = list(duplication = 0.8, death = 0.1, E1 = 0.3,E3=0.6),
  E3 = list(duplication = 0.8, death = 0.1, E1 = 0.3,E2 = 0.2)
)


# add a mutant "A" and set its species rates
rates_genetic_clone_A <- list(
  E1 = list(duplication = 1, death = 0.1,E2 = 0.2,E3 = 0.6),
  E2 = list(duplication = 1, death = 0.1, E1 = 0.3,E3=0.6),
  E3 = list(duplication = 1, death = 0.1, E1 = 0.3,E2=0.2)
)


rates_genetic_clone_B <- list(
  E1 = list(duplication = 1.5, death = 0.1,E2 = 0.2,E3 = 0.6),
  E2 = list(duplication = 1.5, death = 0.1, E1 = 0.3,E3=0.6),
  E3 = list(duplication = 1.5, death = 0.1, E1 = 0.3,E2=0.2)
)

q_G1=get_epigenetic_Q(epigenetic_rates = rates_genetic_clone_0,epistate_name = "G1")
q_G2=get_epigenetic_Q(epigenetic_rates = rates_genetic_clone_A,epistate_name = "G2")
q_G3=get_epigenetic_Q(epigenetic_rates = rates_genetic_clone_B,epistate_name = "G3")

transition_matrix_plot <-q_G1$plot_q
sim$add_mutant("G1", rates_genetic_clone_0)
sim$place_cell("G1[E1]", 500, 500)

# let the simulation evolve until the species "A[E2]" has less than 10 cells
sim$run_up_to_time(20)


sim$add_mutant("G2",rates_genetic_clone_A)
starting_cell=sim$choose_border_cell_in("G1")
plot_tissue(sim)+
  geom_point(data = starting_cell,aes(x=position_x,y=position_y),inherit.aes = FALSE)+
  facet_wrap(~epistate)

sim$mutate_progeny(starting_cell,"G2")
sim$set_rates(list("G1" = list(duplication = 0.6)))
sim$run_up_to_time(40)
plot_tissue(sim)+
  facet_wrap(~epistate)



# Sampling ncells with random box sampling of boxes of size n_w x n_h


sim$add_mutant("G3",rates_genetic_clone_B)
starting_cell=sim$choose_border_cell_in("G2")
plot_tissue(sim)+
  geom_point(data = starting_cell,aes(x=position_x,y=position_y),inherit.aes = FALSE)

sim$mutate_progeny(starting_cell,"G3")
sim$set_rates(list("G2" = list(duplication = 0.6, death = 0.15)))
sim$set_rates(list("G1" = list(duplication = 0.4, death = 0.2)))
sim$run_up_to_time(80)
plot_tissue(sim)+
  facet_wrap(~epistate)

bbox_width=40
bbox1_p <- c(350,430)
bbox1_q <- bbox1_p + bbox_width
# view the boxes
plot_tissue(sim) +
  geom_rect(
    xmin = bbox1_p[1],
    xmax = bbox1_q[1],
    ymin = bbox1_p[2],
    ymax = bbox1_q[2],
    fill = NA,
    color = "black"
  )


sim$sample_cells("S3", bbox1_p, bbox1_q)
sample_forest <- sim$get_sample_forest()


sampled_epigenetic_prop= sample_forest$get_nodes() %>%
  filter(!is.na(sample)) %>%
  group_by(epistate,mutant) %>%
  summarise(n=n())
sampled_epigenetic_prop <- sampled_epigenetic_prop %>%
  group_by(mutant) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

sankey_epi <-ggplot(sampled_epigenetic_prop,
       aes(x = mutant, stratum = epistate, alluvium = epistate,
           y = prop, fill = epistate)) +
  geom_flow(alpha = 0.4) +
  geom_stratum(width = 1/3, color = "grey30") +
  geom_text(stat = "stratum", aes(label = epistate), size = 3) +
  labs(title = "Epistate proportions across mutant status",
       x = "Mutant", y = "Count", fill = "Epistate") +
  scale_fill_manual(values=c("E1"="forestgreen","E2"="goldenrod","E3"="orchid2"))+
  theme_minimal()


sampled_genetic_prop= sample_forest$get_nodes() %>%
  filter(!is.na(sample)) %>%
  group_by(epistate,mutant) %>%
  summarise(n=n())
sampled_genetic_prop <- sampled_genetic_prop %>%
  group_by(epistate) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

sankey_gen <-ggplot(sampled_genetic_prop,
       aes(x = epistate, stratum = mutant, alluvium = mutant,
           y = prop, fill = mutant)) +
  geom_flow(alpha = 0.4) +
  geom_stratum(width = 1/3, color = "grey30") +
  geom_text(stat = "stratum", aes(label = mutant), size = 3) +
  labs(title = "Mutant proportions across epistate",
       x = "Epistate", y = "Count", fill = "Mutant") +
  scale_fill_manual(values=c("G1"="coral2","G2"="turquoise4","G3"="darkorange"))+
  theme_minimal()
library(grid)
ht <- grid.grabExpr(
  draw(transition_matrix_plot)
)
wrap_plots(
  list(ht, sankey_gen, sankey_epi),
  design = "AABB\nCC##"
) + plot_annotation(tag_levels = "a")&
  theme(legend.position = "bottom")
  

color_map_epi = c("G1[E1]"="forestgreen",
                  "G2[E1]"="forestgreen",
                  "G3[E1]"="forestgreen",
                  "G1[E2]"="goldenrod",
                  "G2[E2]"="goldenrod",
                  "G3[E2]"="goldenrod",
                  "G1[E3]"="orchid2",
                  "G2[E3]"="orchid2",
                  "G3[E3]"="orchid2"
                  )
color_map_gen = c("G1[E1]"="coral2",
                  "G2[E1]"="turquoise4",
                  "G3[E1]"="darkorange",
                  "G1[E2]"="coral2",
                  "G2[E2]"="turquoise4",
                  "G3[E2]"="darkorange",
                  "G1[E3]"="coral2",
                  "G2[E3]"="turquoise4",
                  "G3[E3]"="darkorange"
)
p1_forest<-plot_forest(sample_forest,color_map = color_map_epi)+theme(legend.position = "none")+ggtitle(label="Cell coloured by epistate")
p2_forest <-plot_forest(sample_forest,color_map = color_map_gen)+theme(legend.position = "none")+ggtitle(label="Cell coloured by genetic clone")
wrap_plots(
  list(p2_forest,sankey_gen, p1_forest,sankey_epi),
  design = "AABB\nCCDD"
) + plot_annotation(tag_levels = "a")&
  theme(legend.position = "none")
sample_forest$save("sample_forest_atac_case_1.sff")
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

