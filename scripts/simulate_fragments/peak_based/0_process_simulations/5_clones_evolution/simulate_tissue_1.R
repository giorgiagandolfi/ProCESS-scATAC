
rm(list = ls())

library(ProCESS)
library(dplyr)
library(ggplot2)
library(patchwork)

source("../../utils.R")

set.seed(123)

# ============================================================
# 1. Output directories
# ============================================================

# dir.create("results", showWarnings = FALSE)
# dir.create(
#   "results/figures",
#   recursive = TRUE,
#   showWarnings = FALSE
# )

# ============================================================
# 2. Simulation settings
# ============================================================

epistates <- c("P1", "P2", "P3")

# P1 = Differentiated intestinal epithelial
# P2 = Stem-like / proliferative
# P3 = Regenerative / invasive

final_time <- 120

sim <- TissueSimulation(
  epigenetic_states = epistates,
  width = 1000,
  height = 1000,
  save_directory = TRUE,
  name = "CRC_five_clones"
)

# ============================================================
# 3. Helper: clone-specific rates
# ============================================================

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

# ============================================================
# 4. Define the five genetic clones
# ============================================================

# G1: APC loss
# Ancestral clone

rates_G1 <- make_rates()

# G2: APC + KRAS G12D
# Increased P2 proliferation

rates_G2 <- make_rates(
  division = c(
    P1 = 0.8,
    P2 = 1.3,
    P3 = 0.7
  )
)

# G3: APC + KRAS + TP53 loss
# Switching-driven P3 enrichment

rates_G3 <- make_rates(
  division = c(
    P1 = 0.8,
    P2 = 0.9,
    P3 = 0.7
  ),
  p23 = 0.25
)

# G4: APC + KRAS + SMAD4 loss
# Proliferation-driven P3 enrichment

rates_G4 <- make_rates(
  division = c(
    P1 = 0.8,
    P2 = 0.9,
    P3 = 1.4
  )
)

# G5: APC + KRAS + PIK3CA H1047R
# Combined switching and P3 proliferation

rates_G5 <- make_rates(
  division = c(
    P1 = 0.8,
    P2 = 0.9,
    P3 = 1.1
  ),
  p23 = 0.12
)

clone_rates <- list(
  G1 = rates_G1,
  G2 = rates_G2,
  G3 = rates_G3,
  G4 = rates_G4,
  G5 = rates_G5
)

# ============================================================
# 5. Plot epigenetic transition matrices
# ============================================================

q_plots <- lapply(
  names(clone_rates),
  function(g) {
    get_epigenetic_Q(
      epigenetic_rates = clone_rates[[g]],
      epistate_name = g
    )$plot_q
  }
)

names(q_plots) <- names(clone_rates)

p_Q <- wrap_plots(q_plots, ncol = 3)

# ggsave(
#   "results/figures/transition_matrices.pdf",
#   plot = p_Q,
#   width = 13,
#   height = 8
# )

# ============================================================
# 6. Population recording helpers
# ============================================================

population_history <- list()

record_population <- function(sim, label) {
  
  x <- sim$get_counts() %>%
    mutate(
      time = sim$get_clock(),
      stage = label
    )
  
  population_history[[length(population_history) + 1]] <<- x
  
  clone_summary <- x %>%
    group_by(mutant) %>%
    summarise(
      n_cells = sum(counts),
      .groups = "drop"
    ) %>%
    mutate(
      fraction = n_cells / sum(n_cells)
    )
  
  cat(
    "\n========== ",
    label,
    " | time = ",
    round(sim$get_clock(), 2),
    " ==========\n",
    sep = ""
  )
  
  print(clone_summary)
  
  invisible(clone_summary)
}

# Helper: inspect founder availability

check_epistate <- function(sim, clone, epistate) {
  
  n <- sim$get_counts() %>%
    filter(
      mutant == clone,
      .data$epistate == .env$epistate
    ) %>%
    summarise(n = sum(counts)) %>%
    pull(n)
  
  if (length(n) == 0 || is.na(n) || n == 0) {
    stop(
      paste0(
        "No cells available for ",
        clone, "[", epistate, "]"
      )
    )
  }
  
  invisible(n)
}

# ============================================================
# 7. G1: ancestral APC clone
# ============================================================

sim$add_mutant("G1", rates_G1)

sim$place_cell(
  "G1[P1]",
  500,
  500
)

sim$run_up_to_time(15)

record_population(sim, "Before_G2")

# ============================================================
# 8. G2: KRAS clone at t = 15
# ============================================================

sim$add_mutant("G2", rates_G2)

founder <- sim$choose_border_cell_in("G1")

sim$mutate_progeny(
  founder,
  "G2"
)

record_population(sim, "After_G2")

sim$run_up_to_time(30)

record_population(sim, "Before_G3")

# ============================================================
# 9. Reduce ancestral clone growth at t = 30
# ============================================================

# These are synthetic interventions intended
# to prevent G1 and G2 from dominating.

sim$set_rates(list(
  G1 = list(
    duplication = 0.3,
    death = 0.15
  )
))

sim$set_rates(list(
  G2 = list(
    duplication = 0.5,
    death = 0.15
  )
))

# ============================================================
# 10. G3: TP53 branch at t = 30
# ============================================================

sim$add_mutant("G3", rates_G3)

founder <- sim$choose_border_cell_in("G2")

sim$mutate_progeny(
  founder,
  "G3"
)

record_population(sim, "After_G3")

sim$run_up_to_time(35)

record_population(sim, "Before_G4")

# ============================================================
# 11. G4: SMAD4 branch at t = 35
# ============================================================

# G4 benefits from high P3 duplication.
# Start from a G2[P3] founder.

sim$add_mutant("G4", rates_G4)

check_epistate(sim, "G2", "P3")

founder <- sim$choose_border_cell_in("G2[P3]")

if (is.null(founder)) {
  stop("No eligible G2[P3] border founder for G4")
}

sim$mutate_progeny(
  founder,
  "G4"
)

record_population(sim, "After_G4")

sim$run_up_to_time(40)

record_population(sim, "Before_G5")

# ============================================================
# 12. G5: PIK3CA branch at t = 40
# ============================================================

# G5 benefits from both P2 -> P3 switching
# and increased P3 duplication.
#
# Start from a G2[P3] founder to improve
# the chance of early establishment.

sim$add_mutant("G5", rates_G5)

check_epistate(sim, "G2", "P3")

founder <- sim$choose_border_cell_in("G2[P3]")

if (is.null(founder)) {
  stop("No eligible G2[P3] border founder for G5")
}

sim$mutate_progeny(
  founder,
  "G5"
)

record_population(sim, "After_G5")

# ============================================================
# 13. Continue evolution
# ============================================================

# Keep G2 duplication at 0.5.
# Do NOT reduce it further to 0.3.

# Reduce G1 slightly more after all branches exist.

sim$set_rates(list(
  G1 = list(
    duplication = 0.2,
    death = 0.15
  )
))

# ============================================================
# 14. Track population evolution
# ============================================================

time_points <- c(
  45, 50, 60, 70, 80,
  90, 100, 110, 120
)

for (t in time_points) {
  
  sim$run_up_to_time(t)
  
  record_population(
    sim,
    paste0("Time_", t)
  )
}

counts <- sim$get_counts()

phenotype_summary <- counts %>%
  group_by(mutant) %>%
  mutate(
    clone_total = sum(counts),
    fraction = if_else(
      clone_total > 0,
      counts / clone_total,
      0
    )
  ) %>%
  ungroup() %>%
  dplyr::select(mutant, epistate, counts, fraction)

print(phenotype_summary)



# sampling ncells with random box sampling of boxes of size n_w x n_h
# a bounding box 50x50 with at least 100 cells of species B
n_w <- n_h <- 50
ncells <- 0.8 * n_w * n_h
bbox1 <- sim$search_sample(c("G3" = ncells/4, "G5"=ncells/4), n_w, n_h)
bbox2 <- sim$search_sample(c("G1" = ncells/4,"G2"=ncells/4), n_w, n_h)
bbox3 <- sim$search_sample(c("G4" = ncells/4), n_w, n_h)

sim$sample_cells("S1", bbox1$lower_corner, bbox1$upper_corner)
sim$sample_cells("S2", bbox2$lower_corner, bbox2$upper_corner)
sim$sample_cells("S3", bbox3$lower_corner, bbox3$upper_corner)
sample_forest <- sim$get_sample_forest()
plot_forest(sample_forest)

sim$get_samples_info() %>% pull(tumour_cells) %>% sum()
sample_forest$save("sample_forest.sff")


############## MUTATION ENGINE #################
#dir.create("../process_references_v1.3.5/")
setwd("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project/ProCESS-scATAC/scripts/simulate_fragments/peak_based/0_process_simulations/process_references_v1.3.5/")

m_engine <- MutationEngine(setup_code = "GRCh38",tumour_type = "COADREAD", context_sampling = 20,
                           germline_subject = "NA20514",
                           COSMIC_account = list("email"="giorgia.gandolfi@phd.units.it","password"="2*db!XQ4sgQ!dbg"))


mu_SNV = 1e-9
mu_CNA = 0
mu_INDELs = 1e-9

CNA_Clone_G1 <- ProCESS::CNA(
  type = "D",
  "5",
  from = 107000000,
  len = 1e7
)

# G3: TP53 loss — chromosome 17p13.1
# TP53 is located around chr17:7.67–7.69 Mb

CNA_Clone_G3 <- ProCESS::CNA(
  type = "D",
  "17",
  from = 3000000,
  len = 1e7
)

# G4: SMAD4 loss — chromosome 18q21.2
# SMAD4 is located around chr18:51.0–51.1 Mb

CNA_Clone_G4 <- ProCESS::CNA(
  type = "D",
  "18",
  from = 46000000,
  len = 1e7
)

## MYC AMP
CNA_CloneG3_1 = ProCESS::CNA(type = "A", "8",
                            from = 123000000, len = 1e7,src_allele = 1)
CNA_CloneG3_2 = ProCESS::CNA(type = "A", "8",
                            from = 123000000, len = 1e7,src_allele = 1)
CNA_CloneG3_3 = ProCESS::CNA(type = "A", "8",
                            from = 123000000, len = 1e7,src_allele = 1)
CNA_CloneG3_4 = ProCESS::CNA(type = "A", "8",
                            from = 123000000, len = 1e7,src_allele = 1)
CNA_CloneG3_5 = ProCESS::CNA(type = "A", "8",
                            from = 123000000, len = 1e7,src_allele = 1)
CNA_CloneG3_6 = ProCESS::CNA(type = "A", "8",
                            from = 123000000, len = 1e7,src_allele = 1)
CNA_CloneG3_7 = ProCESS::CNA(type = "A", "8",
                            from = 123000000, len = 1e7,src_allele = 1)
CNA_CloneG3_8 = ProCESS::CNA(type = "A", "8",
                            from = 123000000, len = 1e7,src_allele = 1)


## Drivers for the tumors
m_engine$add_mutant(mutant_name = "G1",
                    passenger_rates = list("P1" = c(SNV = mu_SNV, indel = mu_INDELs,CNA=mu_CNA),
                                           "P2" = c(SNV = mu_SNV, indel = mu_INDELs,CNA=mu_CNA),
                                           "P3" = c(SNV = mu_SNV, indel = mu_INDELs,CNA=mu_CNA)),
                    drivers = list(CNA_Clone_G1))
m_engine$add_mutant(mutant_name = "G2",
                    passenger_rates = list("P1" = c(SNV = mu_SNV, indel = mu_INDELs,CNA=mu_CNA),
                                           "P2" = c(SNV = mu_SNV, indel = mu_INDELs,CNA=mu_CNA),
                                           "P3" = c(SNV = mu_SNV, indel = mu_INDELs,CNA=mu_CNA)),
                    drivers = list("KRAS G12D"))

m_engine$add_mutant(mutant_name = "G3",
                    passenger_rates = list("P1" = c(SNV = mu_SNV, indel = mu_INDELs,CNA=mu_CNA),
                                           "P2" = c(SNV = mu_SNV, indel = mu_INDELs,CNA=mu_CNA),
                                           "P3" = c(SNV = mu_SNV, indel = mu_INDELs,CNA=mu_CNA)),
                    drivers = list(CNA_Clone_G3,CNA_CloneG3_1,CNA_CloneG3_2,CNA_CloneG3_3,CNA_CloneG3_4,
                                   CNA_CloneG3_5,CNA_CloneG3_6,CNA_CloneG3_7,CNA_CloneG3_8))

m_engine$add_mutant(mutant_name = "G4",
                    passenger_rates = list("P1" = c(SNV = mu_SNV, indel = mu_INDELs,CNA=mu_CNA),
                                           "P2" = c(SNV = mu_SNV, indel = mu_INDELs,CNA=mu_CNA),
                                           "P3" = c(SNV = mu_SNV, indel = mu_INDELs,CNA=mu_CNA)),
                    drivers = list(CNA_Clone_G4))

m_engine$add_mutant(mutant_name = "G5",
                    passenger_rates = list("P1" = c(SNV = mu_SNV, indel = mu_INDELs,CNA=mu_CNA),
                                           "P2" = c(SNV = mu_SNV, indel = mu_INDELs,CNA=mu_CNA),
                                           "P3" = c(SNV = mu_SNV, indel = mu_INDELs,CNA=mu_CNA)
                                           ),
                    drivers = list("PIK3CA H1047R"))


m_engine$add_exposure(time = 0,coefficients = c(SBS1 = 0.15,SBS5 = 0.40,
                                                SBS18 = 0.15,SBS17b = 0.20,ID1 = 0.40,ID2 = 0.40,ID18=0.2,SBS88 = 0.10))
phylo_forest <- m_engine$place_mutations(sample_forest, num_of_preneoplatic_SNVs=800, num_of_preneoplatic_indels=200)
plot_forest(phylo_forest) %>%
  annotate_forest(phylo_forest,drivers = T,add_driver_label = T)


phylo_forest$save("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project/ProCESS-scATAC/scripts/simulate_fragments/peak_based/0_process_simulations/5_clones_evolution/phylo_forest.sff")


