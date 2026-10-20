# =============================================================================
# CASE IV — REPEATED TREATMENT WITH / WITHOUT EPIGENETIC MEMORY
#
# States:
#   S = sensitive
#   R = resistant / tolerant
#
# Design:
#   Baseline -> Treatment 1 -> Holiday -> Treatment 2
#
# MEMORY CASE:
#   T1 S -> R = 0.30
#   T2 S -> R = 2.00
#
# NO-MEMORY CONTROL:
#   T1 S -> R = 0.30
#   T2 S -> R = 0.30
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
# 1. MODEL INFORMATION
# =============================================================================

# -------------------------------------------------------------------------
# MEMORY CASE
# -------------------------------------------------------------------------

MODEL <- "Epigenetic memory"


SIM_NAME <- "repeated_treatment_with_epigenetic_memory"

T2_SWITCHING_RATE <- 2.00


# -------------------------------------------------------------------------
# FOR THE NO-MEMORY SCRIPT USE INSTEAD:
#
# MODEL <- "No memory"
#
# SIM_NAME <- "repeated_treatment_without_epigenetic_memory"
#
# T2_SWITCHING_RATE <- 0.30
# -------------------------------------------------------------------------


# =============================================================================
# 2. COLOURS
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
# 3. RATES
# =============================================================================


# -----------------------------------------------------------------------------
# Baseline
# -----------------------------------------------------------------------------

# rates_baseline <- list(
#   
#   S = list(
#     duplication = 1.0,
#     death       = 0.6,
#     R           = 0.00001
#   ),
#   
#   R = list(
#     duplication = 1.0,
#     death       = 0.6,
#     S           = 0.10
#   )
# )
# 
# 
# # -----------------------------------------------------------------------------
# # Treatment 1
# # -----------------------------------------------------------------------------
# 
# rates_T1 <- list(
#   
#   S = list(
#     duplication = 1.0,
#     death       = 1.5,
#     R           = 0.30
#   ),
#   
#   R = list(
#     duplication = 1.0,
#     death       = 0.6,
#     S           = 0.00001
#   )
# )
# 
# 
# # -----------------------------------------------------------------------------
# # Holiday
# # -----------------------------------------------------------------------------
# 
# rates_holiday <- list(
#   
#   S = list(
#     duplication = 1.0,
#     death       = 0.6,
#     R           = 0.00001
#   ),
#   
#   R = list(
#     duplication = 1.0,
#     death       = 0.6,
#     S           = 0.30
#   )
# )
# 
# 
# # -----------------------------------------------------------------------------
# # Treatment 2
# #
# # MEMORY:
# #     S -> R = 2.00
# #
# # NO MEMORY:
# #     S -> R = 0.30
# # -----------------------------------------------------------------------------
# 
# rates_T2 <- list(
#   
#   S = list(
#     duplication = 1.0,
#     death       = 1.5,
#     R           = T2_SWITCHING_RATE
#   ),
#   
#   R = list(
#     duplication = 1.0,
#     death       = 0.6,
#     S           = 0.00001
#   )
# )
# 

rates_baseline <- list(
  
  S = list(
    duplication = 1.0,
    death       = 0.75,
    R           = 0.00001
  ),
  
  R = list(
    duplication = 1.0,
    death       = 0.75,
    S           = 0.10 #0.3
  )
)

rates_T1 <- list(
  
  S = list(
    duplication = 1.0,
    death       = 1.30,
    R           = 0.30
  ),
  
  R = list(
    duplication = 1.0,
    death       = 0.75,
    S           = 0.00001
  )
)


rates_holiday <- list(
  
  S = list(
    duplication = 1.0,
    death       = 0.75,
    R           = 0.00001
  ),
  
  R = list(
    duplication = 1.0,
    death       = 0.75,
    S           = 1.00
  )
)


rates_T2 <- list(
  
  S = list(
    duplication = 1.0,
    death       = 1.30,
    
    # Memory effect
    R           = 1.00
  ),
  
  R = list(
    duplication = 1.0,
    death       = 0.75,
    S           = 0.00001
  )
)
# =============================================================================
# 4. PARAMETER SUMMARY
# =============================================================================

parameter_summary <- tibble::tribble(
  
  ~phase,        ~S_duplication, ~S_death, ~S_to_R,
  ~R_duplication, ~R_death, ~R_to_S,
  
  "Baseline",        1.0, 0.6, 0.00001,
  1.0, 0.6, 0.10,
  
  "Treatment 1",     1.0, 1.5, 0.30,
  1.0, 0.6, 0.00001,
  
  "Holiday",         1.0, 0.6, 0.00001,
  1.0, 0.6, 0.30,
  
  "Treatment 2",     1.0, 1.5, T2_SWITCHING_RATE,
  1.0, 0.6, 0.00001
) %>%
  mutate(
    model = MODEL,
    .before = 1
  )



# =============================================================================
# 5. INITIALISE SIMULATION
# =============================================================================

sim <- TissueSimulation(
  epigenetic_states = c("S", "R"),
  width = 2000,
  height = 2000,
  save_directory = TRUE,
  name = SIM_NAME
)

sim$death_activation_level <- 50


sim$add_mutant(
  "G1",
  rates_baseline
)


sim$place_cell(
  "G1[S]",
  1000,
  1000
)


# =============================================================================
# 6. FUNCTION TO SAVE A SNAPSHOT
# =============================================================================

get_snapshot <- function(
    sim,
    model,
    phase,
    stage,
    treatment = NA_integer_
) {
  
  counts <- sim$get_counts() %>%
    filter(
      mutant == "G1",
      epistate %in% c("S", "R")
    ) %>%
    select(
      epistate,
      counts
    ) %>%
    tidyr::pivot_wider(
      names_from = epistate,
      values_from = counts,
      values_fill = 0
    )
  
  # Make sure columns exist
  if (!"S" %in% names(counts)) counts$S <- 0
  if (!"R" %in% names(counts)) counts$R <- 0
  
  counts %>%
    mutate(
      model = model,
      phase = phase,
      stage = stage,
      treatment = treatment,
      
      time = as.numeric(sim$get_clock()),
      
      total = S + R,
      
      prop_S = S / total,
      prop_R = R / total,
      
      percent_S = 100 * prop_S,
      percent_R = 100 * prop_R
    ) %>%
    select(
      model,
      phase,
      stage,
      treatment,
      time,
      S,
      R,
      total,
      prop_S,
      prop_R,
      percent_S,
      percent_R
    )
}


# =============================================================================
# 7. INITIAL STATE
# =============================================================================

snap_initial <- get_snapshot(
  sim = sim,
  model = MODEL,
  phase = "Baseline",
  stage = "Initial"
)


# =============================================================================
# 8. BASELINE
# =============================================================================

v_s <- sim$var("G1[S]")
v_r <- sim$var("G1[R]")
baseline_stopping_criterion <-
  (v_s + v_r) == 2000


sim$run_until(
  baseline_stopping_criterion
)





# end_size <- 2000
# phase <- "Baseline"
# condition_fun <- condition_total
# treatment_number = NA_integer_
# step = 250

# condition_total <- function(sim, target_size) {
#   
#   phenotypes <- sim$get_epigenetic_states() %>%
#     dplyr::pull(epistate)
#   
#   genotypes <- sim$get_mutants() %>%
#     dplyr::pull(mutant)
#   
#   va <- list()
#   
#   for (g in genotypes) {
#     for (p in phenotypes) {
#       
#       name <- paste0(g, "[", p, "]")
#       
#       va[[name]] <- sim$var(name)
#     }
#   }
#   
#   Reduce(`+`, va) == target_size
# }
# counts_baseline <- run_and_record_size(sim = sim,end_size = 2000,phase = "Baseline",condition_fun = condition_total,step=250)  
# plot_timeseries(sim,color_map = STATE_COLORS)




snap_baseline_end <- get_snapshot(
  sim = sim,
  model = MODEL,
  phase = "Baseline",
  stage = "Baseline end"
)

# =============================================================================
# 9. TREATMENT 1
# =============================================================================

sim$set_rates(
  list(
    G1 = rates_T1
  )
)


time_T1_start <- as.numeric(
  sim$get_clock()
)


snap_T1_start <- get_snapshot(
  sim = sim,
  model = MODEL,
  phase = "Treatment 1",
  stage = "T1 start",
  treatment = 1
)


v_s <- sim$var("G1[S]")
v_r <- sim$var("G1[R]")


# -----------------------------------------------------------------------------
# Stop when R >= 4*S
#
# Approximately:
# R >= 80%
# -----------------------------------------------------------------------------

T1_stopping_criterion <-
  v_r >= (
    v_s +
      v_s +
      v_s +
      v_s
  )


sim$run_until(
  T1_stopping_criterion
)


time_T1_end <- as.numeric(
  sim$get_clock()
)


time_to_relapse_T1 <-
  time_T1_end -
  time_T1_start

snap_T1_end <- get_snapshot(
  sim = sim,
  model = MODEL,
  phase = "Treatment 1",
  stage = "T1 end",
  treatment = 1
)
plot_timeseries(sim,color_map = STATE_COLORS)



# =============================================================================
# 10. HOLIDAY
# =============================================================================

sim$set_rates(
  list(
    G1 = rates_holiday
  )
)


time_holiday_start <- as.numeric(
  sim$get_clock()
)


snap_holiday_start <- get_snapshot(
  sim = sim,
  model = MODEL,
  phase = "Holiday",
  stage = "Holiday start"
)


v_s <- sim$var("G1[S]")
v_r <- sim$var("G1[R]")


# -----------------------------------------------------------------------------
# End holiday when S >= 4*R
#
# Approximately:
# S >= 80%
# -----------------------------------------------------------------------------

holiday_stopping_criterion <-
  v_s >= (
    v_r +
      v_r +
      v_r +
      v_r
  )
holiday_stopping_criterion <- v_r<=snap_holiday_start$S

sim$run_until(
  holiday_stopping_criterion
)


time_holiday_end <- as.numeric(
  sim$get_clock()
)


holiday_duration <-
  time_holiday_end -
  time_holiday_start


snap_holiday_end <- get_snapshot(
  sim = sim,
  model = MODEL,
  phase = "Holiday",
  stage = "Holiday end"
)


plot_timeseries(sim,color_map = STATE_COLORS)+
  geom_vline(xintercept = time_T1_start)+
  geom_vline(xintercept = time_T1_end)


print(
  sim$get_counts()
)


cat(
  "\nHoliday duration:",
  holiday_duration,
  "\n"
)


# =============================================================================
# 11. TREATMENT 2
# =============================================================================

sim$set_rates(
  list(
    G1 = rates_T2
  )
)


time_T2_start <- as.numeric(
  sim$get_clock()
)


snap_T2_start <- get_snapshot(
  sim = sim,
  model = MODEL,
  phase = "Treatment 2",
  stage = "T2 start",
  treatment = 2
)


v_s <- sim$var("G1[S]")
v_r <- sim$var("G1[R]")


# -----------------------------------------------------------------------------
# Stop at the same endpoint as T1:
#
# R >= 4*S
# -----------------------------------------------------------------------------

T2_stopping_criterion <-
  v_r >= (
    v_s +
      v_s +
      v_s +
      v_s
  )

# T2_stopping_criterion <- v_r == snap_T1_end$R
sim$run_until(
  T2_stopping_criterion
)


time_T2_end <- as.numeric(
  sim$get_clock()
)


time_to_relapse_T2 <-
  time_T2_end -
  time_T2_start


snap_T2_end <- get_snapshot(
  sim = sim,
  model = MODEL,
  phase = "Treatment 2",
  stage = "T2 end",
  treatment = 2
)


cat("\n--- END TREATMENT 2 ---\n")

print(
  sim$get_counts()
)


cat(
  "\nTime to relapse T2:",
  time_to_relapse_T2,
  "\n"
)


# =============================================================================
# 12. COMBINE POPULATION SNAPSHOTS
# =============================================================================

population_summary <- bind_rows(
  
  snap_initial,
  
  snap_baseline_end,
  
  snap_T1_start,
  
  snap_T1_end,
  
  snap_holiday_start,
  
  snap_holiday_end,
  
  snap_T2_start,
  
  snap_T2_end
  
)


print(
  population_summary
)

plot_timeseries(sim,color_map = STATE_COLORS)

# =============================================================================
# 13. MODEL-LEVEL SUMMARY STATISTICS
# =============================================================================

summary_stats <- tibble(
  
  model = MODEL,
  
  S_to_R_T1 = rates_T1$S$R,
  
  S_to_R_T2 = rates_T2$S$R,
  
  switching_fold_change =
    rates_T2$S$R /
    rates_T1$S$R,
  
  time_T1_start =
    time_T1_start,
  
  time_T1_end =
    time_T1_end,
  
  time_to_relapse_T1 =
    time_to_relapse_T1,
  
  holiday_start =
    time_holiday_start,
  
  holiday_end =
    time_holiday_end,
  
  holiday_duration =
    holiday_duration,
  
  time_T2_start =
    time_T2_start,
  
  time_T2_end =
    time_T2_end,
  
  time_to_relapse_T2 =
    time_to_relapse_T2
  
)


print(
  summary_stats
)


# =============================================================================
# 14. TREATMENT-SPECIFIC SUMMARY
#
# This table is particularly useful for plotting T1 vs T2.
# =============================================================================

treatment_summary <- tibble(
  
  model = c(
    MODEL,
    MODEL
  ),
  
  treatment = c(
    "T1",
    "T2"
  ),
  
  switching_S_to_R = c(
    rates_T1$S$R,
    rates_T2$S$R
  ),
  
  start_time = c(
    time_T1_start,
    time_T2_start
  ),
  
  end_time = c(
    time_T1_end,
    time_T2_end
  ),
  
  time_to_relapse = c(
    time_to_relapse_T1,
    time_to_relapse_T2
  )
  
)


print(
  treatment_summary
)


# =============================================================================
# 15. SAVE EVERYTHING
# =============================================================================

output_dir <- file.path(
  prj_dir,
  "comparison_memory"
)


dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


# Clean model name for filenames
model_file_name <- ifelse(
  MODEL == "Epigenetic memory",
  "memory",
  "no_memory"
)


# -----------------------------------------------------------------------------
# Parameters
# -----------------------------------------------------------------------------

write.csv(
  parameter_summary,
  file.path(
    output_dir,
    paste0(
      model_file_name,
      "_parameters.csv"
    )
  ),
  row.names = FALSE
)


# -----------------------------------------------------------------------------
# Population composition at important points
# -----------------------------------------------------------------------------

write.csv(
  population_summary,
  file.path(
    output_dir,
    paste0(
      model_file_name,
      "_population_summary.csv"
    )
  ),
  row.names = FALSE
)


# -----------------------------------------------------------------------------
# One row per model
# -----------------------------------------------------------------------------

write.csv(
  summary_stats,
  file.path(
    output_dir,
    paste0(
      model_file_name,
      "_summary_stats.csv"
    )
  ),
  row.names = FALSE
)


# -----------------------------------------------------------------------------
# One row per treatment
# -----------------------------------------------------------------------------

write.csv(
  treatment_summary,
  file.path(
    output_dir,
    paste0(
      model_file_name,
      "_treatment_summary.csv"
    )
  ),
  row.names = FALSE
)


# =============================================================================
# 16. FINAL OUTPUT
# =============================================================================

cat("\n")
cat("============================================================\n")
cat(MODEL, "\n")
cat("============================================================\n")

cat(
  "T1 S -> R:",
  rates_T1$S$R,
  "\n"
)

cat(
  "T2 S -> R:",
  rates_T2$S$R,
  "\n"
)

cat(
  "Switching fold-change:",
  rates_T2$S$R / rates_T1$S$R,
  "\n"
)

cat(
  "T1 time to relapse:",
  time_to_relapse_T1,
  "\n"
)

cat(
  "T2 time to relapse:",
  time_to_relapse_T2,
  "\n"
)

cat("============================================================\n")


print(
  sim$get_counts()
)

