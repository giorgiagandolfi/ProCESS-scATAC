
rm(list = ls())

library(ProCESS)
library(dplyr)
library(ggplot2)
library(patchwork)
source("../utils.R")
set.seed(123)

# dir.create("results", showWarnings = FALSE)
# dir.create("results/figures", recursive = TRUE,
#            showWarnings = FALSE)

# --------------------------------------------------
# 1. Define the three epistates
# --------------------------------------------------

epistates <- c("P1", "P2", "P3")

# P1 = Differentiated epithelial
# P2 = Stem-like / proliferative
# P3 = Regenerative / invasive

sim <- TissueSimulation(
  epigenetic_states = epistates,
  width = 1000,
  height = 1000,
  save_directory = FALSE,
  name = "CRC_five_clones"
)

# --------------------------------------------------
# 2. Helper for clone-specific rates
# --------------------------------------------------

make_rates <- function(
    division = c(P1 = 0.8, P2 = 0.9, P3 = 0.7),
    death = c(P1 = 0.1, P2 = 0.1, P3 = 0.1),
    p12 = 0.03,
    p13 = 0.005,
    p21 = 0.04,
    p23 = 0.01,
    p31 = 0.005,
    p32 = 0.04
) {
  
  list(
    P1 = list(
      duplication = unname(division["P1"]),
      death = unname(death["P1"]),
      P2 = p12,
      P3 = p13
    ),
    
    P2 = list(
      duplication = unname(division["P2"]),
      death = unname(death["P2"]),
      P1 = p21,
      P3 = p23
    ),
    
    P3 = list(
      duplication = unname(division["P3"]),
      death = unname(death["P3"]),
      P1 = p31,
      P2 = p32
    )
  )
}

# --------------------------------------------------
# 3. Five genetic clones
# --------------------------------------------------

# G1: APC loss
# Baseline dynamics

rates_G1 <- make_rates()

# G2: APC + KRAS G12D
# Higher P2 proliferation

rates_G2 <- make_rates(
  division = c(P1 = 0.8, P2 = 1.3, P3 = 0.7)
)

# G3: APC + KRAS + TP53 loss
# P3 enrichment driven by switching

rates_G3 <- make_rates(
  division = c(P1 = 0.8, P2 = 0.9, P3 = 0.7),
  p23 = 0.25
)

# G4: APC + KRAS + SMAD4 loss
# P3 enrichment driven by proliferation

rates_G4 <- make_rates(
  division = c(P1 = 0.8, P2 = 0.9, P3 = 1.4)
)

# G5: APC + KRAS + PIK3CA H1047R
# P3 enrichment driven by switching + proliferation

# rates_G5 <- make_rates(
#   division = c(P1 = 0.8, P2 = 0.9, P3 = 1.4),
#   p23 = 0.25
# )
rates_G5 <- make_rates(
  division = c(P1 = 0.8, P2 = 0.9, P3 = 1.1),
  p23 = 0.12
)
clone_rates <- list(
  G1 = rates_G1,
  G2 = rates_G2,
  G3 = rates_G3,
  G4 = rates_G4,
  G5 = rates_G5
)

# --------------------------------------------------
# 4. Generate transition matrix plots
# --------------------------------------------------

q_plots <- lapply(names(clone_rates), function(g) {
  
  get_epigenetic_Q(
    epigenetic_rates = clone_rates[[g]],
    epistate_name = g
  )$plot_q
  
})

names(q_plots) <- names(clone_rates)

p_Q <- wrap_plots(q_plots, ncol = 3)

# ggsave(
#   "results/figures/transition_matrices.pdf",
#   p_Q,
#   width = 13,
#   height = 8
# )

# --------------------------------------------------
# 5. Ancestral clone: APC
# --------------------------------------------------
sim$add_mutant("G1", rates_G1)
sim$place_cell("G1[P1]", 500, 500)

sim$run_up_to_time(15)

# --------------------------------------------------
# G2: KRAS clone
# --------------------------------------------------

sim$add_mutant("G2", rates_G2)

founder <- sim$choose_border_cell_in("G1")
sim$mutate_progeny(founder, "G2")

sim$run_up_to_time(30)

# Reduce ancestral growth before branching
sim$set_rates(list(
  G1 = list(duplication = 0.3, death = 0.15)
))
sim$set_rates(list(
  G2 = list(duplication = 0.5, death = 0.15)
))

# --------------------------------------------------
# G3: TP53 clone
# --------------------------------------------------

sim$add_mutant("G3", rates_G3)

founder <- sim$choose_border_cell_in("G2")
sim$mutate_progeny(founder, "G3")

sim$run_up_to_time(35)

# --------------------------------------------------
# G4: SMAD4 clone
# --------------------------------------------------
# Reduce ancestral growth before branching
sim$set_rates(list(
  G1 = list(duplication = 0.2, death = 0.15)
))
sim$set_rates(list(
  G2 = list(duplication = 0.4, death = 0.15)
))


sim$add_mutant("G4", rates_G4)

founder <- sim$choose_border_cell_in("G2")
sim$mutate_progeny(founder, "G4")

sim$run_up_to_time(40)
plot_state(sim)
# --------------------------------------------------
# G5: PIK3CA clone
# --------------------------------------------------

sim$add_mutant("G5", rates_G5)

founder <- sim$choose_border_cell_in("G2")
sim$mutate_progeny(founder, "G5")

# # Reduce G2 growth after all branches exist
# sim$set_rates(list(
#   G2 = list(duplication = 0.3, death = 0.15)
# ))
# sim$set_rates(list(
#   G1 = list(duplication = 0.2, death = 0.15)
# ))
# 
# --------------------------------------------------
# Continue evolution
# --------------------------------------------------

sim$run_up_to_time(100)
plot_tissue(sim)+facet_wrap(~mutant)
plot_state(sim)