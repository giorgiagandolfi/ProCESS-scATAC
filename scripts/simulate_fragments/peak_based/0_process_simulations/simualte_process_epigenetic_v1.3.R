rm(list=ls())
library(ProCESS)
library(dplyr)





# set the seed of the random number generator for repeatability
set.seed(0)
mutant_cols <- c(
  "A[E1]" = "goldenrod",
  "A[E2]" = "magenta4",
  "A[E3]" = "royalblue3",
  "A[E4]" ="forestgreen"
)
# create a tissue simulation with two epigenetic states
sim <- TissueSimulation(epigenetic_states = c("E1", "E2","E3","E4"),
                        width = 300, height = 300)

# add a mutant "A" and set its species rates
rates <- list(
  E1 = list(duplication = 0.8, death = 0.1, E2 = 0.5, E3 = 0.4, E4 = 0.01),
  E2 = list(duplication = 0.8, death = 0.1, E1 = 0.12, E3 = 0.3, E4 = 0.2),
  E3 = list(duplication = 0.8, death = 0.1, E1 = 0.4, E2 = 0.3, E4 = 0.02),
  E4 = list(duplication = 0.8, death = 0.1, E1 = 0.01, E2 = 0.8, E3 = 0.02)
)
sim$add_mutant("A", rates)

states <- names(rates)


sim$place_cell("A[E1]", 150, 150)

# let the simulation evolve until the species "A[E2]" has less than 10 cells
sim$run_up_to_time(30)

plot_tissue(simulation = sim,color_map = mutant_cols)
plot_state(sim,color_map = mutant_cols)


n_w <- n_h <- 25
ncells <- 0.8 * n_w * n_h

# Sampling ncells with random box sampling of boxes of size n_w x n_h
bbox <- sim$search_sample(c("A" = ncells), n_w, n_h)
sim$sample_cells("S1", bbox$lower_corner, bbox$upper_corner)
sample_forest <- sim$get_sample_forest()

# plot it
plot_forest(sample_forest,color_map = mutant_cols) %>%
  annotate_forest(sample_forest)



sample_forest$save("sample_forest_atac_epigenome_new_version.sff")
setwd("process_references_v1.3")
m_engine <- MutationEngine(setup_code = "GRCh38",tumour_type = "COADREAD", context_sampling = 20,
                           COSMIC_account = list("email"="giorgia.gandolfi@phd.units.it","password"="2*db!XQ4sgQ!dbg"))


mu_SNV = 1e-9
mu_CNA = 5e-12
mu_INDELs = 1e-9


CNA_Clone1 = ProCESS::CNA(type = "A", "12",
                          from = 10000000, len = 4e7)

## Drivers for the tumors
m_engine$add_mutant(mutant_name = "A",
                    passenger_rates = list("E1" = c(SNV = mu_SNV, indel = mu_INDELs, CNA=mu_CNA),
                                           "E2" = c(SNV = mu_SNV, indel = mu_INDELs, CNA=mu_CNA),
                                           "E3" = c(SNV = mu_SNV, indel = mu_INDELs, CNA=mu_CNA),
                                           "E4" = c(SNV = mu_SNV, indel = mu_INDELs, CNA=1e-11)),
                    drivers = list("KRAS G12D"))

m_engine$add_exposure(time = 0,coefficients = c(SBS1 = 0.15,SBS5 = 0.40,
                                                SBS18 = 0.15,SBS17b = 0.20,ID1 = 0.40,ID2 = 0.40,ID18=0.2,SBS88 = 0.10))
phylo_forest <- m_engine$place_mutations(sample_forest, num_of_preneoplatic_SNVs=800, num_of_preneoplatic_indels=200)
phylo_forest$save("../phylo_forest_atac_epigenetic_new_version.sff")
