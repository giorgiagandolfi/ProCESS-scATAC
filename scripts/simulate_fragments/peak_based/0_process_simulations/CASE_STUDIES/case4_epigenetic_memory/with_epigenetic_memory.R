# =============================================================================
# CASE IV — EPIGENETIC MEMORY UNDER REPEATED TREATMENT
#
# States:
#   S = sensitive
#   R = resistant / tolerant
#
# Design:
#   Baseline -> Treatment 1 -> Holiday -> Treatment 2
#
# Memory mechanism:
#   S -> R switching is higher during T2 than during T1.
# =============================================================================


# =============================================================================
# 0. SETUP
# =============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(ProCESS)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

set.seed(0)

prj_dir <- paste0(
  "/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/",
  "scATAC_project/ProCESS-scATAC/scripts/simulate_fragments/",
  "peak_based/0_process_simulations/CASE_STUDIES/",
  "case4_epigenetic_memory/"
)

setwd(prj_dir)

source("utils_epigenetic_memory.R")


# =============================================================================
# 1. COLOURS
# =============================================================================

STATE_COLORS <- c(
  "G1[S]" = "deepskyblue4",
  "G1[R]" = "firebrick2"
)

EPISTATE_COLORS <- c(
  "S" = "deepskyblue4",
  "R" = "firebrick2"
)


# =============================================================================
# 2. MODEL PARAMETERS
# =============================================================================


# -----------------------------------------------------------------------------
# Baseline
#
# S dominates.
# S -> R is essentially absent.
# R -> S favours return to the sensitive state.
# -----------------------------------------------------------------------------

rates_baseline <- list(
  
  S = list(
    duplication = 1.0,
    death       = 0.6,
    R           = 0.00001
  ),
  
  R = list(
    duplication = 1.0,
    death       = 0.6,
    S           = 0.10
  )
)


# -----------------------------------------------------------------------------
# Treatment 1
#
# Drug suppresses S.
# R is favoured under treatment.
# Treatment also induces S -> R switching.
# -----------------------------------------------------------------------------

rates_T1 <- list(
  
  S = list(
    duplication = 1.0,
    death       = 1.5,
    R           = 0.30
  ),
  
  R = list(
    duplication = 1.0,
    death       = 0.6,
    S           = 0.00001
  )
)


# -----------------------------------------------------------------------------
# Holiday
#
# No fitness difference between S and R.
#
# R decreases because it switches back to S,
# NOT because R has a higher death rate.
# -----------------------------------------------------------------------------

rates_holiday <- list(
  
  S = list(
    duplication = 1.0,
    death       = 0.6,
    R           = 0.00001
  ),
  
  R = list(
    duplication = 1.0,
    death       = 0.6,
    S           = 0.30
  )
)


# -----------------------------------------------------------------------------
# Treatment 2
#
# Same fitness landscape as T1.
#
# ONLY important difference:
#
#       T1: S -> R = 0.30
#       T2: S -> R = 0.60
#
# This increased switching represents epigenetic memory.
# -----------------------------------------------------------------------------

rates_T2 <- list(
  
  S = list(
    duplication = 1.0,
    death       = 1.5,
    R           = 2
  ),
  
  R = list(
    duplication = 1.0,
    death       = 0.6,
    S           = 0.00001
  )
)


# =============================================================================
# 3. PARAMETER SUMMARY
# =============================================================================

parameter_summary <- tibble::tribble(
  
  ~phase,       ~S_duplication, ~S_death, ~S_to_R,
  ~R_duplication, ~R_death, ~R_to_S,
  
  "Baseline",       1.0, 0.6, 0.00001,
  1.0, 0.6, 0.10,
  
  "Treatment 1",    1.0, 1.5, 0.30,
  1.0, 0.6, 0.00001,
  
  "Holiday",        1.0, 0.6, 0.00001,
  1.0, 0.6, 0.30,
  
  "Treatment 2",    1.0, 1.5, 0.60,
  1.0, 0.6, 0.00001
)

print(parameter_summary)


# =============================================================================
# 4. INITIALISE SIMULATION
# =============================================================================

sim <- TissueSimulation(
  epigenetic_states = c("S", "R"),
  width = 2000,
  height = 2000,
  save_directory = TRUE,
  name = "epigenetic_memory_double_treatment_new"
)

sim$death_activation_level <- 50

sim$add_mutant(
  "G1",
  rates_baseline
)

# Start from one sensitive cell
sim$place_cell(
  "G1[S]",
  1000,
  1000
)


# =============================================================================
# 5. BASELINE
# =============================================================================

count_initial <- get_state_counts(
  sim,
  phase = "Baseline"
)


# Variables used by ProCESS/C++
v_s <- sim$var("G1[S]")
v_r <- sim$var("G1[R]")


# Grow baseline tumour to 2,000 cells
baseline_stopping_criterion <-
  (v_s + v_r) == 2000


sim$run_until(
  baseline_stopping_criterion
)


count_baseline_end <- get_state_counts(
  sim,
  phase = "Baseline"
)


cat("\n--- END BASELINE ---\n")
print(sim$get_counts())
plot_timeseries(sim,color_map = STATE_COLORS)

# =============================================================================
# 6. TREATMENT 1
# =============================================================================

sim$set_rates(
  list(G1 = rates_T1)
)


count_T1_start <- get_state_counts(
  sim,
  phase = "Treatment 1",
  treatment_number = 1
)


# Re-create variables for clarity
v_s <- sim$var("G1[S]")
v_r <- sim$var("G1[R]")


# -----------------------------------------------------------------------------
# End T1 when R strongly dominates.
#
# R >= 4*S corresponds approximately to:
#
# R >= 80%
# S <= 20%
# -----------------------------------------------------------------------------

T1_stopping_criterion <-
  v_r >= (v_s + v_s + v_s + v_s)


sim$run_until(
  T1_stopping_criterion
)


count_T1_end <- get_state_counts(
  sim,
  phase = "Treatment 1",
  treatment_number = 1
)


cat("\n--- END TREATMENT 1 ---\n")
print(sim$get_counts())
plot_timeseries(sim,color_map = STATE_COLORS)


# =============================================================================
# 7. HOLIDAY
# =============================================================================

sim$set_rates(
  list(G1 = rates_holiday)
)


count_holiday_start <- get_state_counts(
  sim,
  phase = "Holiday",
  treatment_number = NA_integer_
)


v_s <- sim$var("G1[S]")
v_r <- sim$var("G1[R]")


# -----------------------------------------------------------------------------
# During holiday:
#
# R -> S = 0.30
# S -> R ~ 0
#
# S and R have identical duplication/death rates.
#
# Therefore recovery of S is driven by phenotype switching.
#
# Stop when S >= 4*R:
#
# S ~ 80%
# R ~ 20%
# -----------------------------------------------------------------------------

holiday_stopping_criterion <-
  v_s >= (v_r + v_r + v_r + v_r)


sim$run_until(
  holiday_stopping_criterion
)


count_holiday_end <- get_state_counts(
  sim,
  phase = "Holiday",
  treatment_number = NA_integer_
)


cat("\n--- END HOLIDAY ---\n")
print(sim$get_counts())
plot_timeseries(sim,color_map = STATE_COLORS)


# =============================================================================
# 8. TREATMENT 2
# =============================================================================

sim$set_rates(
  list(G1 = rates_T2)
)


count_T2_start <- get_state_counts(
  sim,
  phase = "Treatment 2",
  treatment_number = 2
)


v_s <- sim$var("G1[S]")
v_r <- sim$var("G1[R]")


# -----------------------------------------------------------------------------
# Same endpoint as T1.
#
# Critical difference:
#
# T1: S -> R = 0.30
# T2: S -> R = 0.60
#
# Therefore R should emerge more rapidly during the second exposure.
# -----------------------------------------------------------------------------

T2_stopping_criterion <-
  v_r >= (v_s + v_s + v_s + v_s)


sim$run_until(
  T2_stopping_criterion
)


count_T2_end <- get_state_counts(
  sim,
  phase = "Treatment 2",
  treatment_number = 2
)


cat("\n--- END TREATMENT 2 ---\n")
print(sim$get_counts())
plot_timeseries(sim,color_map = STATE_COLORS)
sim$get_clock()
