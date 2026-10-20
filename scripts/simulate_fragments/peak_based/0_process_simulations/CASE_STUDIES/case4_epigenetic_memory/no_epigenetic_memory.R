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



MODEL <- "No memory"

SIM_NAME <- "repeated_treatment_without_epigenetic_memory"

T2_SWITCHING_RATE <- 0.10

STATE_COLORS <- c(
  "G1[S]" = "deepskyblue4",
  "G1[R]" = "firebrick2"
)

EPISTATE_COLORS <- c(
  "S" = "deepskyblue4",
  "R" = "firebrick2"
)
plot_firings_my <- function(simulation){
  {
    firings <- simulation$get_firings() %>% ProCESS:::add_species_col()
    plot <- ggplot2::ggplot(firings) + ggplot2::geom_bar(stat = "identity", 
                                                         ggplot2::aes(x = "", y = .data$fired, fill = .data$event))
    if ("epistate" %in% colnames(firings)) {
      plot <- plot + ggplot2::facet_grid(.data$mutant ~ .data$epistate)
    }
    else {
      plot <- plot + ggplot2::facet_grid(.data$mutant)
    }
    plot + ggplot2::labs(fill = "Event",x = NULL, y = NULL) + ProCESS:::my_theme() + ggplot2::scale_fill_brewer(palette = "Dark2") + 
      ggplot2::theme(legend.position = "bottom")
  }
}


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

rates_baseline <- list(
  
  S = list(
    duplication = 1.0,
    death       = 0.75,
    R           = 0.00001
  ),
  
  R = list(
    duplication = 1.0,
    death       = 0.75,
    S           = 0.1 #0.30
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
    R           = 0.1
  ),
  
  R = list(
    duplication = 1.0,
    death       = 0.75,
    S           = 0.00001
  )
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


snap_baseline_end <- get_snapshot(
  sim = sim,
  model = MODEL,
  phase = "Baseline",
  stage = "Baseline end"
)

p_state_1 <- plot_state(sim, color_map = STATE_COLORS)+ggtitle(label = "Baseline")
p_firings_1 <- plot_firings_my(sim)+ggtitle(label = "Baseline")

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

p_state_2 <- plot_state(sim, color_map = STATE_COLORS)+ggtitle(label = "Treatment 1")
p_firings_2 <- plot_firings_my(sim)+ggtitle(label = "Treatment 1")

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

cat(
  "\nHoliday duration:",
  holiday_duration,
  "\n"
)

p_state_3 <- plot_state(sim, color_map = STATE_COLORS)+ggtitle(label = "Holiday")
p_firings_3 <- plot_firings_my(sim)+ggtitle(label = "Holiday")


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
p_state_4 <- plot_state(sim, color_map = STATE_COLORS)+ggtitle(label = "Treatment 2")
p_firings_4 <- plot_firings_my(sim)+ggtitle(label = "Treatment 2")





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
# =============================================================================
# EXTRACT TREATMENT TIMES
# =============================================================================

T1_start <- population_summary %>%
  filter(stage == "T1 start") %>%
  pull(time)

T1_end <- population_summary %>%
  filter(stage == "T1 end") %>%
  pull(time)

T2_start <- population_summary %>%
  filter(stage == "T2 start") %>%
  pull(time)

T2_end <- population_summary %>%
  filter(stage == "T2 end") %>%
  pull(time)
p_timeseries <- plot_timeseries(sim, color_map = STATE_COLORS)+
  annotate(
    "rect",
    xmin = T1_start,
    xmax = T1_end,
    ymin = -Inf,
    ymax = Inf,
    fill = "grey50",
    alpha = 0.15
  ) +
  
  annotate(
    "rect",
    xmin = T2_start,
    xmax = T2_end,
    ymin = -Inf,
    ymax = Inf,
    fill = "grey50",
    alpha = 0.15
  ) 

wrap_plots(list(p_state_1,p_state_2,p_state_3,p_state_4,p_firings_1,p_firings_2,p_firings_3,p_firings_4,p_timeseries),
           design = "ABCD\nEFGH\nEFGH\nIIII",guides = "collect") & theme(legend.position = "bottom")

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

