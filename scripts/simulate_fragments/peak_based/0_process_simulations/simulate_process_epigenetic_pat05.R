rm(list=ls())
library(ProCESS)
library(dplyr)
library(ggplot2)
library(tidyverse)



# set the seed of the random number generator for repeatability
set.seed(0)
mutant_cols <- c(
  "A[E1]" = "goldenrod",
  "A[E2]" = "magenta4",
  "A[E3]" = "royalblue3",
  "A[E4]" ="forestgreen",
  "B[E1]" = "goldenrod",
  "B[E2]" = "magenta4",
  "B[E3]" = "royalblue3",
  "B[E4]" ="forestgreen"
)
# create a tissue simulation with two epigenetic states
sim <- TissueSimulation(epigenetic_states = c("E1", "E2","E3","E4"),
                        width = 1000, height = 1000,save_snapshots = F)

rates_genetic_clone_0 <- list(
  E1 = list(duplication = 0.6, death = 0.1),
  E2 = list(duplication = 0.6, death = 0.1),
  E3 = list(duplication = 0.6, death = 0.1),
  E4 = list(duplication = 0.6, death = 0.1)
)


# add a mutant "A" and set its species rates
rates_genetic_clone_A <- list(
  E1 = list(duplication = 0.8, death = 0.1, E2 = 0.3, E3 = 0.5, E4 = 0.01),
  E2 = list(duplication = 0.8, death = 0.1, E1 = 0.12, E3 = 0.3, E4 = 0.2),
  E3 = list(duplication = 0.8, death = 0.1, E1 = 0.01, E2 = 0.2, E4 = 0.02),
  E4 = list(duplication = 0.8, death = 0.1, E1 = 0.02, E2 = 0.2, E3 = 0.5)
)

rates_genetic_clone_B <- list(
  E1 = list(duplication = 1, death = 0.1, E2 = 0.3, E3 = 0.5, E4 = 0.1),
  E2 = list(duplication = 0.8, death = 0.1, E1 = 0.8, E3 = 0.5, E4 = 0.2),
  E3 = list(duplication = 0.8, death = 0.1, E1 = 0.4, E2 = 0.2, E4 = 0.1),
  E4 = list(duplication = 0.8, death = 0.1, E1 = 0.2, E2 = 0.5, E3 = 0.02)
)

rates_genetic_clone_C <- list(
  E1 = list(duplication = 1.5, death = 0.1, E2 = 0.3, E3 = 0.2, E4 = 0.1),
  E2 = list(duplication = 1.5, death = 0.1, E1 = 0.8, E3 = 0.1, E4 = 0.1),
  E3 = list(duplication = 1.4, death = 0.1, E1 = 0.4, E2 = 0.2, E4 = 0.1),
  E4 = list(duplication = 1.4, death = 0.1, E1 = 0.4, E2 = 0.5, E3 = 0.02)
)

rates_genetic_clone_D <- list(
  E1 = list(duplication = 2, death = 0.1, E2 = 0.3, E3 = 0.2, E4 = 0.1),
  E2 = list(duplication = 2, death = 0.1, E1 = 0.8, E3 = 0.1, E4 = 0.1),
  E3 = list(duplication = 1.6, death = 0.1, E1 = 0.4, E2 = 0.2, E4 = 0.1),
  E4 = list(duplication = 1.6, death = 0.1, E1 = 0.4, E2 = 0.5, E3 = 0.02)
)


sim$add_mutant("0", rates_genetic_clone_0)

# states <- names(rates_genetic_clone_B)
#
# # Initialize matrix
# mat <- matrix(0,
#               nrow = length(states),
#               ncol = length(states),
#               dimnames = list(states, states))
#
# # Fill in switching rates
# for (from in states) {
#   for (to in names(rates_genetic_clone_B[[from]])) {
#     if (!(to %in% c("duplication", "death"))) {
#       mat[from, to] <- rates_genetic_clone_B[[from]][[to]]
#     }
#   }
# }
# diag(mat) <- -rowSums(mat - diag(diag(mat)))
#
#
#
# Heatmap(mat,col = colorRampPalette(c("white", "steelblue"))(10),cluster_rows = F,cluster_columns = F,
#         cell_fun = function(j, i, x, y, width, height, fill) {
#           grid.text(sprintf("%.2f", mat[i, j]), x, y, gp = gpar(fontsize = 10))
#         })
# place one cell of "A" in epigenetic state "E1"
sim$place_cell("0[E1]", 500, 500)

# let the simulation evolve until the species "A[E2]" has less than 10 cells
sim$run_up_to_time(20)
plot_tissue(simulation = sim)
plot_state(sim)


sim$add_mutant("A",rates_genetic_clone_A)
# choose a border cell in "A" and let one of its progeny mutate in "B"
# sim$mutate_progeny(sim$choose_border_cell_in("A"), "B")
sim$mutate_progeny(sim$choose_border_cell_in("0[E1]"),"A")
sim$set_rates(list("0" = list(duplication = 0.4)))
sim$run_up_to_time(40)
plot_state(sim)

n_w <- n_h <- 15
ncells <- 0.8 * n_w * n_h

# Sampling ncells with random box sampling of boxes of size n_w x n_h
bbox <- sim$search_sample(c("A" = ncells), n_w, n_h)
sim$sample_cells("S1", bbox$lower_corner, bbox$upper_corner)

sim$add_mutant("B",rates_genetic_clone_B)
sim$mutate_progeny(sim$choose_border_cell_in("0[E1]"),"B")
sim$set_rates(list("A" = list(duplication = 0.4, death = 0.1)))
sim$set_rates(list("0" = list(duplication = 0.05, death = 0.8)))
sim$run_up_to_time(60)
plot_state(sim)



sim$add_mutant("C",rates_genetic_clone_C)
sim$mutate_progeny(sim$choose_border_cell_in("B[E1]"),"C")
sim$set_rates(list("A" = list(duplication = 0.3, death = 0.8)))
sim$set_rates(list("B" = list(duplication = 0.4, death = 0.8)))
sim$run_up_to_time(80)
plot_state(sim)


sim$add_mutant("D",rates_genetic_clone_D)
sim$mutate_progeny(sim$choose_border_cell_in("C[E1]"),"D")
sim$set_rates(list("A" = list(duplication = 0.2, death = .9)))
sim$set_rates(list("B" = list(duplication = 0.2, death = .9)))
sim$set_rates(list("C" = list(duplication = 0.4, death = .2)))
sim$run_up_to_time(100)
plot_state(sim)
plot_tissue(sim)


n_w <- n_h <- 15
ncells <- 0.8 * n_w * n_h

# Sampling ncells with random box sampling of boxes of size n_w x n_h
bbox <- sim$search_sample(c("C" = ncells/2,"D"=ncells/2), n_w, n_h)
sim$sample_cells("S2", bbox$lower_corner, bbox$upper_corner)

bbox <- sim$search_sample(c("B" = ncells), n_w, n_h)
sim$sample_cells("S3", bbox$lower_corner, bbox$upper_corner)

sample_forest <- sim$get_sample_forest()

# plot it
plot_forest(sample_forest) %>%
  annotate_forest(sample_forest)



sample_forest$save("sample_forest_atac_epigenome_new_version_pat05.sff")
setwd("process_references_v1.3")

m_engine <- MutationEngine(setup_code = "GRCh38",tumour_type = "COADREAD", context_sampling = 20,
                           germline_subject = "NA20514",
                           COSMIC_account = list("email"="giorgia.gandolfi@phd.units.it","password"="2*db!XQ4sgQ!dbg"))


mu_SNV = 1e-9
mu_CNA = 5e-12
mu_INDELs = 1e-9


# pat_05_cn = readRDS("final_pat05_copy_number.rds")

CNA_Clone0_1 = ProCESS::CNA(type = "D", "11",
                            from = 48500001, len = 1e7)
CNA_Clone0_2  = ProCESS::CNA(type = "D", "5",
                             from = 1000001, len = 109499999,allele = 0)
CNA_Clone0_3 = ProCESS::CNA(type = "D","17",
                            from = 500001, len = 42999999)

# chrX      5 107000001 150500000       0.25         1 chrX:107000001:150500000 43499999 D

CNA_CloneA = ProCESS::CNA(type = "D", "X",
                          from = 107000001, len = 43499999)



CNA_CloneB = ProCESS::CNA(type = "A", "8",
                          from = 127118340, len = 1e7)


CNA_CloneC = ProCESS::CNA(type = "D", "8",
                          from = 107000001, len = 22999999)

CNA_CloneD = ProCESS::CNA(type = "D", "14",
                          from = 46500001, len = 22999999)

## Drivers for the tumors
m_engine$add_mutant(mutant_name = "0",
                    passenger_rates = list("E1" = c(SNV = mu_SNV, indel = mu_INDELs, CNA=0),
                                           "E2" = c(SNV = mu_SNV, indel = mu_INDELs, CNA=0),
                                           "E3" = c(SNV = mu_SNV, indel = mu_INDELs, CNA=0),
                                           "E4" = c(SNV = mu_SNV, indel = mu_INDELs, CNA=0)),
                    drivers = list(list("APC R1450*", allele = 1),CNA_Clone0_2,CNA_Clone0_3))
m_engine$add_mutant(mutant_name = "A",
                    passenger_rates = list("E1" = c(SNV = mu_SNV, indel = mu_INDELs, CNA=0),
                                           "E2" = c(SNV = mu_SNV, indel = mu_INDELs, CNA=0),
                                           "E3" = c(SNV = mu_SNV, indel = mu_INDELs, CNA=0),
                                           "E4" = c(SNV = mu_SNV, indel = mu_INDELs, CNA=0)),
                    drivers = list(CNA_CloneA))

m_engine$add_mutant(mutant_name = "B",
                    passenger_rates = list("E1" = c(SNV = mu_SNV, indel = mu_INDELs, CNA=0),
                                           "E2" = c(SNV = mu_SNV, indel = mu_INDELs, CNA=0),
                                           "E3" = c(SNV = mu_SNV, indel = mu_INDELs, CNA=0),
                                           "E4" = c(SNV = mu_SNV, indel = mu_INDELs, CNA=0)),
                    drivers = list(CNA_CloneB))

m_engine$add_mutant(mutant_name = "C",
                    passenger_rates = list("E1" = c(SNV = mu_SNV, indel = mu_INDELs, CNA=0),
                                           "E2" = c(SNV = mu_SNV, indel = mu_INDELs, CNA=0),
                                           "E3" = c(SNV = mu_SNV, indel = mu_INDELs, CNA=0),
                                           "E4" = c(SNV = mu_SNV, indel = mu_INDELs, CNA=0)),
                    drivers = list(CNA_CloneC))


m_engine$add_mutant(mutant_name = "D",
                    passenger_rates = list("E1" = c(SNV = mu_SNV, indel = mu_INDELs, CNA=0),
                                           "E2" = c(SNV = mu_SNV, indel = mu_INDELs, CNA=0),
                                           "E3" = c(SNV = mu_SNV, indel = mu_INDELs, CNA=0),
                                           "E4" = c(SNV = mu_SNV, indel = mu_INDELs, CNA=0)),
                    drivers = list(CNA_CloneD))
m_engine$add_exposure(time = 0,coefficients = c(SBS1 = 0.15,SBS5 = 0.40,
                                                SBS18 = 0.15,SBS17b = 0.20,ID1 = 0.40,ID2 = 0.40,ID18=0.2,SBS88 = 0.10))
phylo_forest <- m_engine$place_mutations(sample_forest, num_of_preneoplatic_SNVs=800, num_of_preneoplatic_indels=200)
phylo_forest$save("../phylo_forest_atac_epigenetic_new_version_pat05.sff")

