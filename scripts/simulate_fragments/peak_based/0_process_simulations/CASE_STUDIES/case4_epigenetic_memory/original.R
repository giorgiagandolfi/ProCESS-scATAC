rm(list = ls())

library(ProCESS)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggalluvial)
library(patchwork)
library(scales)

source(
  "/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project/ProCESS-scATAC/scripts/simulate_fragments/peak_based/utils.R"
)

set.seed(0)


############################################################
# EPIGENETIC MEMORY MODEL
#
# P1 = drug-sensitive phenotype
# P2 = drug-resistant / drug-tolerant phenotype
#
#
# TIMELINE
#
# 0   -------- 80   Baseline
# 80  -------- 120  Treatment 1
# 120 -------- 160  Recovery / relapse
# 160 -------- 200  Treatment 2
# 200 -------- 240  Recovery / relapse
#
#
# T1 duration = T2 duration = 40
#
#
# KEY HYPOTHESIS
#
# Drug pressure is IDENTICAL during T1 and T2.
#
# Duplication/death rates:
#       T1 == T2
#
# But:
#
# P1 -> P2 switching
#
# T1 = 0.10
# T2 = 0.30
#
# The increased switching during T2 represents
# EPIGENETIC MEMORY.
#
############################################################



############################################################
# 1. GLOBAL TIME PARAMETERS
############################################################

baseline_end <- 80

treatment_duration <- 40
recovery_duration <- 40


# Treatment 1

T1_start <- baseline_end
T1_end <- T1_start + treatment_duration


# Recovery 1

recovery1_start <- T1_end
recovery1_end <- recovery1_start + recovery_duration


# Treatment 2

T2_start <- recovery1_end
T2_end <- T2_start + treatment_duration


# Recovery 2

recovery2_start <- T2_end
recovery2_end <- recovery2_start + recovery_duration


############################################################
# Record every 2 simulation units
#
# This gives better temporal resolution of killing
# and switching dynamics.
############################################################

record_every <- 2



############################################################
# 2. INITIALIZE SIMULATION
############################################################

sim <- TissueSimulation(
  epigenetic_states = c("P1", "P2"),
  width = 2000,
  height = 2000,
  save_directory = TRUE,
  name = "epigenetic_memory_double_treatment"
)

sim$death_activation_level <- 50



############################################################
# 3. COLORS
############################################################

state_colors <- c(
  "G1[P1]" = "goldenrod",
  "G1[P2]" = "orchid2"
)

epistate_colors <- c(
  "P1" = "goldenrod",
  "P2" = "orchid2"
)



############################################################
# 4. FUNCTION TO RECORD CELL COUNTS
############################################################

get_state_counts <- function(
    sim,
    phase,
    treatment_number = NA
) {
  
  sim$get_cells() %>%
    
    dplyr::group_by(epistate) %>%
    
    dplyr::summarise(
      n = n(),
      .groups = "drop"
    ) %>%
    
    dplyr::mutate(
      
      total_cells = sum(n),
      
      prop = n / total_cells,
      
      time = sim$get_clock(),
      
      phase = phase,
      
      treatment = treatment_number
    )
}



############################################################
# 5. FUNCTION TO RUN AND RECORD SIMULATION
############################################################

run_and_record <- function(
    sim,
    end_time,
    phase,
    treatment_number = NA,
    step = 2
) {
  
  results <- list()
  
  current_time <- sim$get_clock()
  
  recording_times <- seq(
    from = current_time + step,
    to = end_time,
    by = step
  )
  
  
  # Make sure final time is included
  
  if (length(recording_times) == 0) {
    
    recording_times <- end_time
    
  } else if (tail(recording_times, 1) < end_time) {
    
    recording_times <- c(
      recording_times,
      end_time
    )
  }
  
  
  for (t in recording_times) {
    
    sim$run_up_to_time(t)
    
    results[[length(results) + 1]] <-
      
      get_state_counts(
        sim = sim,
        phase = phase,
        treatment_number = treatment_number
      )
  }
  
  
  bind_rows(results)
}



############################################################
# 6. FUNCTION TO SAMPLE TISSUE
############################################################

sample_tissue <- function(
    sim,
    sample_name,
    n_w = 25,
    n_h = 25
) {
  
  ncells <- 0.8 * n_w * n_h
  
  
  bbox <- sim$search_sample(
    c("G1" = ncells),
    n_w,
    n_h
  )
  
  
  sim$sample_cells(
    sample_name,
    bbox$lower_corner,
    bbox$upper_corner
  )
  
  
  p <- plot_tissue(
    sim,
    at_sample = sample_name,
    color_map = state_colors
  ) +
    
    facet_wrap(~epistate) +
    
    ggtitle(sample_name)
  
  
  return(p)
}



############################################################
# 7. BASELINE RATES
#
# P1:
# dominant phenotype without treatment.
#
# P2:
# resistant phenotype with some fitness cost.
#
#
# P1 -> P2 = rare spontaneous switching
#
# P2 -> P1 = relatively frequent reversion
#
############################################################

rates_baseline <- list(
  
  P1 = list(
    
    duplication = 1.0,
    
    death = 0.10,
    
    # P1 -> P2
    P2 = 0.01
  ),
  
  
  P2 = list(
    
    duplication = 0.80,
    
    death = 0.10,
    
    # P2 -> P1
    P1 = 0.30
  )
)



############################################################
# 8. TREATMENT 1
#
# Mild / gradual killing.
#
#
# P1:
#
# duplication = 0.45
# death       = 0.47
#
# P1 is only mildly negatively selected.
#
#
# P2:
#
# duplication = 0.45
# death       = 0.30
#
# P2 survives better under treatment.
#
#
# P1 -> P2 = 0.10
#
############################################################

rates_T1 <- list(
  
  P1 = list(
    
    duplication = 0.45,
    
    death = 0.47,
    
    # Sensitive -> resistant
    P2 = 0.10
  ),
  
  
  P2 = list(
    
    duplication = 0.45,
    
    death = 0.30,
    
    # Resistant -> sensitive
    P1 = 0.02
  )
)



############################################################
# 9. RECOVERY RATES
#
# Drug is removed.
#
# P1 regains its proliferative advantage.
#
# P2 can revert toward P1.
#
# Increasing P2 -> P1 during recovery helps restore
# a predominantly P1 population before T2.
#
############################################################

rates_recovery <- list(
  
  P1 = list(
    
    duplication = 1.00,
    
    death = 0.10,
    
    # Low spontaneous acquisition of P2
    P2 = 0.02
  ),
  
  
  P2 = list(
    
    duplication = 0.80,
    
    death = 0.10,
    
    # Reversion toward sensitive phenotype
    P1 = 0.40
  )
)



############################################################
# 10. TREATMENT 2
#
# IMPORTANT:
#
# DUPLICATION AND DEATH ARE IDENTICAL TO T1.
#
#
# P1:
#
# T1 duplication = 0.45
# T2 duplication = 0.45
#
# T1 death = 0.47
# T2 death = 0.47
#
#
# P2:
#
# T1 duplication = 0.45
# T2 duplication = 0.45
#
# T1 death = 0.30
# T2 death = 0.30
#
#
# ONLY MAJOR DIFFERENCE:
#
# P1 -> P2
#
# T1 = 0.10
# T2 = 0.30
#
# This represents epigenetic memory.
#
############################################################

rates_T2 <- list(
  
  P1 = list(
    
    # SAME AS T1
    duplication = 0.45,
    
    # SAME AS T1
    death = 0.47,
    
    # EPIGENETIC MEMORY
    P2 = 0.30
  ),
  
  
  P2 = list(
    
    # SAME AS T1
    duplication = 0.45,
    
    # SAME AS T1
    death = 0.30,
    
    # SAME AS T1
    P1 = 0.02
  )
)



############################################################
# 11. ADD GENETIC CLONE
############################################################

sim$add_mutant(
  "G1",
  rates_baseline
)



############################################################
# 12. PLACE INITIAL CELL
############################################################

sim$place_cell(
  "G1[P1]",
  1000,
  1000
)



############################################################
# 13. RECORD INITIAL CONDITION
############################################################

count_initial <- get_state_counts(
  sim = sim,
  phase = "Baseline"
)



############################################################
# 14. BASELINE EXPANSION
#
# TIME:
#
# 0 -> 80
#
############################################################

count_baseline <- run_and_record(
  
  sim = sim,
  
  end_time = baseline_end,
  
  phase = "Baseline",
  
  step = record_every
)



############################################################
# PRE-TREATMENT STATE
############################################################

p_pre <- plot_state(
  sim,
  color_map = state_colors
) +
  
  ggtitle(
    "Pre-treatment"
  )



############################################################
# SAMPLE PRE-T1
############################################################

sample_pre <- sample_tissue(
  sim,
  "S_PRE_T1"
)



############################################################
# 15. TREATMENT 1
#
# TIME:
#
# 80 -> 120
#
# Duration = 40
#
############################################################

sim$set_rates(
  list(
    "G1" = rates_T1
  )
)



############################################################
# Record T1 starting condition
############################################################

count_T1_start <- get_state_counts(
  
  sim,
  
  phase = "Treatment 1",
  
  treatment_number = 1
)



############################################################
# Run T1
############################################################

count_T1 <- run_and_record(
  
  sim = sim,
  
  end_time = T1_end,
  
  phase = "Treatment 1",
  
  treatment_number = 1,
  
  step = record_every
)



############################################################
# POST-T1 STATE
############################################################

p_T1 <- plot_state(
  sim,
  color_map = state_colors
) +
  
  ggtitle(
    "Post-Treatment 1"
  )



############################################################
# SAMPLE POST-T1
############################################################

sample_post_T1 <- sample_tissue(
  sim,
  "S_POST_T1"
)



############################################################
# 16. RECOVERY 1
#
# TIME:
#
# 120 -> 160
#
# Duration = 40
#
############################################################

sim$set_rates(
  list(
    "G1" = rates_recovery
  )
)



############################################################
# Run recovery
############################################################

count_recovery1 <- run_and_record(
  
  sim = sim,
  
  end_time = recovery1_end,
  
  phase = "Recovery 1",
  
  step = record_every
)



############################################################
# PRE-T2 STATE
############################################################

p_pre_T2 <- plot_state(
  sim,
  color_map = state_colors
) +
  
  ggtitle(
    "Pre-Treatment 2 / Relapse"
  )



############################################################
# SAMPLE PRE-T2
############################################################

sample_pre_T2 <- sample_tissue(
  sim,
  "S_PRE_T2"
)



############################################################
# 17. TREATMENT 2
#
# TIME:
#
# 160 -> 200
#
# Duration = 40
#
#
# Same treatment pressure as T1.
#
# Higher P1 -> P2 switching represents memory.
#
############################################################

sim$set_rates(
  list(
    "G1" = rates_T2
  )
)



############################################################
# Record T2 starting condition
############################################################

count_T2_start <- get_state_counts(
  
  sim,
  
  phase = "Treatment 2",
  
  treatment_number = 2
)



############################################################
# Run T2
############################################################

count_T2 <- run_and_record(
  
  sim = sim,
  
  end_time = T2_end,
  
  phase = "Treatment 2",
  
  treatment_number = 2,
  
  step = record_every
)



############################################################
# POST-T2 STATE
############################################################

p_T2 <- plot_state(
  sim,
  color_map = state_colors
) +
  
  ggtitle(
    "Post-Treatment 2"
  )



############################################################
# SAMPLE POST-T2
############################################################

sample_post_T2 <- sample_tissue(
  sim,
  "S_POST_T2"
)



############################################################
# 18. RECOVERY 2
#
# TIME:
#
# 200 -> 240
#
# Same recovery parameters as after T1.
#
############################################################

sim$set_rates(
  list(
    "G1" = rates_recovery
  )
)



############################################################
# Run second recovery
############################################################

count_recovery2 <- run_and_record(
  
  sim = sim,
  
  end_time = recovery2_end,
  
  phase = "Recovery 2",
  
  step = record_every
)



############################################################
# FINAL STATE
############################################################

p_final <- plot_state(
  sim,
  color_map = state_colors
) +
  
  ggtitle(
    "Second relapse"
  )



############################################################
# FINAL SAMPLE
############################################################

sample_final <- sample_tissue(
  sim,
  "S_RELAPSE_T2"
)



############################################################
# 19. FINAL COUNTS
############################################################

sim$get_counts()



############################################################
# 20. COMBINE ALL TEMPORAL DATA
############################################################

count_time <- bind_rows(
  
  count_initial,
  
  count_baseline,
  
  count_T1_start,
  
  count_T1,
  
  count_recovery1,
  
  count_T2_start,
  
  count_T2,
  
  count_recovery2
)



############################################################
# Convert time to numeric
############################################################

count_time <- count_time %>%
  
  mutate(
    
    time = as.numeric(
      as.character(time)
    ),
    
    epistate = factor(
      epistate,
      levels = c(
        "P1",
        "P2"
      )
    )
  )


print(count_time)



############################################################
# 21. PARAMETER SUMMARY
############################################################

parameter_summary <- tibble(
  
  phase = c(
    "Baseline",
    "Treatment 1",
    "Recovery 1",
    "Treatment 2",
    "Recovery 2"
  ),
  
  
  P1_duplication = c(
    1.00,
    0.45,
    1.00,
    0.45,
    1.00
  ),
  
  
  P1_death = c(
    0.10,
    0.47,
    0.10,
    0.47,
    0.10
  ),
  
  
  P1_to_P2 = c(
    0.01,
    0.10,
    0.02,
    0.30,
    0.02
  ),
  
  
  P2_duplication = c(
    0.80,
    0.45,
    0.80,
    0.45,
    0.80
  ),
  
  
  P2_death = c(
    0.10,
    0.30,
    0.10,
    0.30,
    0.10
  ),
  
  
  P2_to_P1 = c(
    0.30,
    0.02,
    0.40,
    0.02,
    0.40
  )
)


print(parameter_summary)



############################################################
# 22. EXTRACT TOTAL POPULATION
#
# Important:
#
# We need absolute cell numbers to visualize treatment
# killing.
#
############################################################

total_population <- count_time %>%
  
  select(
    time,
    phase,
    total_cells
  ) %>%
  
  distinct() %>%
  
  arrange(time)



############################################################
# 23. TOTAL CELL NUMBER THROUGH TIME
############################################################

p_population <- ggplot(
  total_population,
  aes(
    x = time,
    y = total_cells
  )
) +
  
  ##########################################################
# Treatment 1 window
##########################################################

annotate(
  "rect",
  xmin = T1_start,
  xmax = T1_end,
  ymin = -Inf,
  ymax = Inf,
  alpha = 0.10
) +
  
  
  ##########################################################
# Treatment 2 window
##########################################################

annotate(
  "rect",
  xmin = T2_start,
  xmax = T2_end,
  ymin = -Inf,
  ymax = Inf,
  alpha = 0.10
) +
  
  
  geom_line(
    linewidth = 1.2
  ) +
  
  
  geom_point(
    size = 1.8
  ) +
  
  
  geom_vline(
    xintercept = c(
      T1_start,
      T1_end,
      T2_start,
      T2_end
    ),
    linetype = "dashed",
    linewidth = 0.4
  ) +
  
  
  annotate(
    "text",
    x = mean(c(T1_start, T1_end)),
    y = Inf,
    label = "Treatment 1",
    vjust = 1.5,
    fontface = "bold"
  ) +
  
  
  annotate(
    "text",
    x = mean(c(T2_start, T2_end)),
    y = Inf,
    label = "Treatment 2",
    vjust = 1.5,
    fontface = "bold"
  ) +
  
  
  labs(
    x = "Simulation time",
    y = "Total number of cells",
    title = "Tumor burden during repeated treatment"
  ) +
  
  
  theme_classic(
    base_size = 14
  )


p_population



############################################################
# 24. ABSOLUTE P1/P2 CELL NUMBERS
#
# This is probably the most informative plot for checking
# whether treatment killing is sufficiently gradual.
############################################################

p_absolute <- plot_timeseries(sim,color_map = state_colors) +
  
  ##########################################################
# Treatment windows
##########################################################

annotate(
  "rect",
  xmin = T1_start,
  xmax = T1_end,
  ymin = -Inf,
  ymax = Inf,
  alpha = 0.10
) +
  
  annotate(
    "rect",
    xmin = T2_start,
    xmax = T2_end,
    ymin = -Inf,
    ymax = Inf,
    alpha = 0.10
  ) +

  
  labs(
    x = "Simulation time",
    y = "Cell count",
    color = "Phenotype",
    title = "Sensitive and resistant populations",
    subtitle = "Mild treatment pressure produces gradual population changes"
  ) +
  
  
  theme_minimal(
    base_size = 14
  ) +
  
  
  theme(
    legend.position = "top"
  )


p_absolute



############################################################
# 25. COMPOSITION / ALLUVIAL-LIKE AREA PLOT
############################################################

p_composition <- ggplot(
  count_time,
  aes(
    x = time,
    y = prop,
    fill = epistate
  )
) +
  
  ##########################################################
# Treatment windows
##########################################################

annotate(
  "rect",
  xmin = T1_start,
  xmax = T1_end,
  ymin = -Inf,
  ymax = Inf,
  alpha = 0.10
) +
  
  annotate(
    "rect",
    xmin = T2_start,
    xmax = T2_end,
    ymin = -Inf,
    ymax = Inf,
    alpha = 0.10
  ) +
  
  
  geom_area(
    position = "stack",
    alpha = 0.85,
    linewidth = 0.3
  ) +
  
  
  scale_fill_manual(
    values = epistate_colors,
    labels = c(
      "P1" = "Sensitive (P1)",
      "P2" = "Resistant (P2)"
    )
  ) +
  
  
  scale_y_continuous(
    labels = scales::percent_format(),
    limits = c(0, 1),
    expand = c(0, 0)
  ) +
  
  
  geom_vline(
    xintercept = c(
      T1_start,
      T1_end,
      T2_start,
      T2_end
    ),
    linetype = "dashed",
    linewidth = 0.4
  ) +
  
  
  labs(
    x = "Simulation time",
    y = "Cell population",
    fill = "Phenotype",
    title = "Epigenetic state composition",
    subtitle = "P1 = sensitive; P2 = resistant"
  ) +
  
  
  theme_classic(
    base_size = 14
  ) +
  
  
  theme(
    legend.position = "top"
  )


p_composition



############################################################
# 26. DIRECT T1 vs T2 COMPARISON
#
# Convert absolute simulation time into:
#
# time since treatment started
#
# T1: 0 -> 40
# T2: 0 -> 40
#
############################################################

treatment_dynamics <- bind_rows(
  
  count_T1_start,
  
  count_T1,
  
  count_T2_start,
  
  count_T2
  
) %>%
  
  mutate(
    
    time = as.numeric(
      as.character(time)
    ),
    
    treatment_time = case_when(
      
      treatment == 1 ~ time - T1_start,
      
      treatment == 2 ~ time - T2_start
    ),
    
    treatment_label = case_when(
      
      treatment == 1 ~ "Treatment 1",
      
      treatment == 2 ~ "Treatment 2"
    )
  )



############################################################
# 27. COMPARE P2 ACQUISITION
############################################################

P2_dynamics <- treatment_dynamics %>%
  
  filter(
    epistate == "P2"
  )



p_P2_memory <- ggplot(
  P2_dynamics,
  aes(
    x = treatment_time,
    y = prop,
    linetype = treatment_label,
    group = treatment_label
  )
) +
  
  geom_line(
    linewidth = 1.3
  ) +
  
  geom_point(
    size = 2
  ) +
  
  scale_y_continuous(
    labels = scales::percent_format(),
    limits = c(0, 1)
  ) +
  
  scale_x_continuous(
    breaks = seq(
      0,
      treatment_duration,
      by = 5
    )
  ) +
  
  labs(
    x = "Time since treatment start",
    y = "P2 resistant-cell proportion",
    linetype = "Exposure",
    title = "Acquisition of the resistant phenotype",
    subtitle = "T1: P1→P2 = 0.10; T2: P1→P2 = 0.30"
  ) +
  
  theme_classic(
    base_size = 14
  )


p_P2_memory



############################################################
# 28. COMPARE TOTAL CELL NUMBER DURING T1 vs T2
#
# Useful for determining whether the two treatments have
# similar killing kinetics.
############################################################

treatment_population <- treatment_dynamics %>%
  
  select(
    treatment_time,
    treatment_label,
    total_cells
  ) %>%
  
  distinct()



p_treatment_killing <- ggplot(
  treatment_population,
  aes(
    x = treatment_time,
    y = total_cells,
    linetype = treatment_label
  )
) +
  
  geom_line(
    linewidth = 1.3
  ) +
  
  geom_point(
    size = 2
  ) +
  
  scale_x_continuous(
    breaks = seq(
      0,
      treatment_duration,
      by = 5
    )
  ) +
  
  labs(
    x = "Time since treatment start",
    y = "Total number of cells",
    linetype = "Exposure",
    title = "Treatment response",
    subtitle = "Comparison of tumor burden during first and second exposure"
  ) +
  
  theme_classic(
    base_size = 14
  )


p_treatment_killing



############################################################
# 29. NORMALIZED TUMOR BURDEN
#
# This is better than absolute cell count for comparing
# treatment response because T1 and T2 can start with
# different numbers of cells.
#
# Each treatment starts at 100%.
#
############################################################

normalized_treatment_population <- treatment_population %>%
  
  group_by(
    treatment_label
  ) %>%
  
  arrange(
    treatment_time,
    .by_group = TRUE
  ) %>%
  
  mutate(
    
    initial_cells = first(total_cells),
    
    relative_burden = total_cells / initial_cells
  ) %>%
  
  ungroup()



p_normalized_killing <- ggplot(
  normalized_treatment_population,
  aes(
    x = treatment_time,
    y = relative_burden,
    linetype = treatment_label
  )
) +
  
  geom_hline(
    yintercept = 1,
    linetype = "dotted"
  ) +
  
  geom_line(
    linewidth = 1.3
  ) +
  
  geom_point(
    size = 2
  ) +
  
  scale_y_continuous(
    labels = scales::percent_format()
  ) +
  
  scale_x_continuous(
    breaks = seq(
      0,
      treatment_duration,
      by = 5
    )
  ) +
  
  labs(
    x = "Time since treatment start",
    y = "Relative tumor burden",
    linetype = "Exposure",
    title = "Normalized treatment response",
    subtitle = "Tumor burden at treatment start = 100%"
  ) +
  
  theme_classic(
    base_size = 14
  )


p_normalized_killing



############################################################
# 30. SPATIAL EVOLUTION
############################################################

p_spatial <- (
  
  p_pre
  
) | (
  
  p_T1
  
) | (
  
  p_pre_T2
  
) | (
  
  p_T2
  
) | (
  
  p_final
  
)


p_spatial



############################################################
# 31. SAMPLE COMPARISON
############################################################

sample_plots <- (
  
  sample_pre
  
) | (
  
  sample_post_T1
  
) | (
  
  sample_pre_T2
  
) | (
  
  sample_post_T2
  
)


sample_plots



############################################################
# 32. FINAL SUMMARY FIGURE
#
# Panel 1:
# absolute tumor burden
#
# Panel 2:
# absolute P1/P2 populations
#
# Panel 3:
# relative composition
#
# Panel 4:
# direct P2 memory comparison
#
############################################################

final_plot <- (
  
  # p_population /
    
    p_absolute /
    
    p_composition /
    
    p_P2_memory
  
) +
  
  plot_layout(
    heights = c(
      # 1,
      1,
      1,
      1
    )
  )


final_plot



############################################################
# 33. SECOND SUMMARY:
#
# DIRECT COMPARISON OF T1 AND T2
############################################################

treatment_comparison_plot <- (
  
  p_normalized_killing /
    
    p_P2_memory
  
) +
  
  plot_layout(
    heights = c(
      1,
      1
    )
  )


treatment_comparison_plot

P2_dynamics_normalized <- treatment_dynamics %>%
  
  filter(epistate == "P2") %>%
  
  group_by(treatment_label) %>%
  
  arrange(treatment_time, .by_group = TRUE) %>%
  
  mutate(
    
    P2_start = first(prop),
    
    delta_P2 = prop - P2_start
    
  ) %>%
  
  ungroup()
p_memory_delta <- ggplot(
  P2_dynamics_normalized,
  aes(
    x = treatment_time,
    y = delta_P2,
    linetype = treatment_label
  )
) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dotted"
  ) +
  
  geom_line(
    linewidth = 1.3
  ) +
  
  geom_point(
    size = 2
  ) +
  
  scale_y_continuous(
    labels = scales::percent_format()
  ) +
  
  labs(
    x = "Time since treatment start",
    y = "Change in P2 proportion",
    linetype = "Treatment",
    title = "P2 acquisition during repeated treatment",
    subtitle = "Each treatment is normalized to its starting P2 proportion"
  ) +
  
  theme_classic(base_size = 14)

p_memory_delta



library(dplyr)
library(ggplot2)
library(patchwork)

############################################################
# 34. PLOT RATE HISTORY
############################################################

rates_history <- sim$get_rates_update_history()

rates_history <- rates_history %>%
  mutate(
    time = as.numeric(time),
    epistate = factor(
      epistate,
      levels = c("P1", "P2")
    )
  )

print(rates_history)


############################################################
# Treatment intervals
############################################################

T1_start <- 80
T1_end   <- 120

T2_start <- 160
T2_end   <- 200


############################################################
# Common phenotype colors
############################################################

epistate_colors <- c(
  "P1" = "goldenrod",
  "P2" = "orchid2"
)


############################################################
# 2. DUPLICATION + DEATH RATES
############################################################

rates_growth <- rates_history %>%
  filter(
    event %in% c("duplication", "death")
  )


p_rates_growth <- ggplot(
  rates_growth,
  aes(
    x = time,
    y = rate,
    color = epistate,
    linetype = event,
    group = interaction(epistate, event)
  )
) +
  
  ##########################################################
# Treatment windows
##########################################################

annotate(
  "rect",
  xmin = T1_start,
  xmax = T1_end,
  ymin = -Inf,
  ymax = Inf,
  # fill = "firebrick",
  alpha = 0.08
) +
  
  annotate(
    "rect",
    xmin = T2_start,
    xmax = T2_end,
    ymin = -Inf,
    ymax = Inf,
    # fill = "forestgreen",
    alpha = 0.08
  ) +
  
  ##########################################################
# Rates
##########################################################

geom_step(
  linewidth = 1.2,
  direction = "hv"
) +
  
  geom_point(
    size = 2.5
  ) +
  
  ##########################################################
# Treatment boundaries
##########################################################

geom_vline(
  xintercept = c(
    T1_start,
    T1_end,
    T2_start,
    T2_end
  ),
  linetype = "dotted",
  linewidth = 0.4
) +
  
  ##########################################################
# Colors
##########################################################

scale_color_manual(
  values = epistate_colors,
  labels = c(
    "P1" = "P1 sensitive",
    "P2" = "P2 resistant"
  )
) +
  
  ##########################################################
# Axes
##########################################################

scale_x_continuous(
  breaks = seq(0, 240, 20),
  limits = c(0, 240),
  expand = c(0, 0)
) +
  
  ##########################################################
# Labels
##########################################################

labs(
  x = "Simulation time",
  y = "Event rate",
  color = "Phenotype",
  linetype = "Event",
  title = "Duplication and death rates",
  subtitle = "Treatment 1 and Treatment 2 have identical fitness parameters"
) +
  
  theme_minimal(
    base_size = 14
  ) +
  
  theme(
    legend.position = "top",
    plot.title = element_text(
      face = "bold"
    )
  )


p_rates_growth


############################################################
# 3. SWITCHING RATES
############################################################

rates_switch <- rates_history %>%
  filter(
    event == "switch"
  ) %>%
  
  mutate(
    transition = case_when(
      epistate == "P1" ~ "P1 → P2",
      epistate == "P2" ~ "P2 → P1"
    )
  )


############################################################
# Plot
############################################################

p_rates_switch <- ggplot(
  rates_switch,
  aes(
    x = time,
    y = rate,
    color = transition,
    group = transition
  )
) +
  
  ##########################################################
# Treatment 1
##########################################################

annotate(
  "rect",
  xmin = T1_start,
  xmax = T1_end,
  ymin = -Inf,
  ymax = Inf,
  # fill = "firebrick",
  alpha = 0.08
) +
  
  ##########################################################
# Treatment 2
##########################################################

annotate(
  "rect",
  xmin = T2_start,
  xmax = T2_end,
  ymin = -Inf,
  ymax = Inf,
  # fill = "forestgreen",
  alpha = 0.08
) +
  
  ##########################################################
# Switching rates
##########################################################

geom_step(
  linewidth = 1.3,
  direction = "hv"
) +
  
  geom_point(
    size = 2.7
  ) +
  
  ##########################################################
# Treatment boundaries
##########################################################

geom_vline(
  xintercept = c(
    T1_start,
    T1_end,
    T2_start,
    T2_end
  ),
  linetype = "dotted",
  linewidth = 0.4
) +
  
  ##########################################################
# Axis
##########################################################

scale_x_continuous(
  breaks = seq(0, 240, 20),
  limits = c(0, 240),
  expand = c(0, 0)
) +
  
  scale_y_continuous(
    limits = c(0, 0.45),
    breaks = seq(0, 0.4, 0.1)
  ) +
  
  ##########################################################
# Labels
##########################################################

labs(
  x = "Simulation time",
  y = "Switching rate",
  color = "Transition",
  title = "Epigenetic switching rates",
  subtitle = "Higher P1 → P2 switching during the second treatment represents memory"
) +
  
  theme_minimal(
    base_size = 14
  ) +
  
  theme(
    legend.position = "top",
    plot.title = element_text(
      face = "bold"
    )
  )


p_rates_switch

