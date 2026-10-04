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
setwd("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project/ProCESS-scATAC/scripts/simulate_fragments/peak_based/0_process_simulations/process_epigenetic_memory")
prj_dir <- "/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project/ProCESS-scATAC/scripts/simulate_fragments/peak_based/0_process_simulations/process_epigenetic_memory"

source("utils_epigenetic_memory.R")
set.seed(0)

# -----------------------------------------------------------------------------
# 1. Experimental design
# -----------------------------------------------------------------------------
# P1 = drug-sensitive phenotype
# P2 = drug-resistant / drug-tolerant phenotype
#
# T1 and T2 have the same duration and the same duplication/death rates.
# The memory hypothesis is encoded only by a higher P1 -> P2 switching rate in T2.

treatment_duration <- 40
recovery_duration <- 40
record_every <- 2

TIMELINE <- list(
  baseline_end = 80,
  T1_start = 80,
  T1_end = 80 + treatment_duration,
  recovery1_end = 80 + treatment_duration + recovery_duration,
  T2_start = 80 + treatment_duration + recovery_duration,
  T2_end = 80 + 2 * treatment_duration + recovery_duration,
  recovery2_end = 80 + 2 * treatment_duration + 2 * recovery_duration
)

STATE_COLORS <- c("G1[P1]" = "goldenrod", "G1[P2]" = "orchid2")
EPISTATE_COLORS <- c("P1" = "goldenrod", "P2" = "orchid2")

# -----------------------------------------------------------------------------
# 2. Rates
# -----------------------------------------------------------------------------
rates_baseline <- list(
  P1 = list(duplication = 1.00, death = 0.10, P2 = 0.01),
  P2 = list(duplication = 0.80, death = 0.10, P1 = 0.30)
)

rates_T1 <- list(
  P1 = list(duplication = 0.45, death = 0.47, P2 = 0.10),
  P2 = list(duplication = 0.45, death = 0.30, P1 = 0.02)
)

rates_recovery <- list(
  P1 = list(duplication = 1.00, death = 0.10, P2 = 0.02),
  P2 = list(duplication = 0.80, death = 0.10, P1 = 0.40)
)

rates_T2 <- list(
  P1 = list(duplication = 0.45, death = 0.47, P2 = 0.30),
  P2 = list(duplication = 0.45, death = 0.30, P1 = 0.02)
)

parameter_summary <- tibble::tribble(
  ~phase,        ~P1_duplication, ~P1_death, ~P1_to_P2, ~P2_duplication, ~P2_death, ~P2_to_P1,
  "Baseline",              1.00,      0.10,       0.01,            0.80,      0.10,       0.30,
  "Treatment 1",           0.45,      0.47,       0.10,            0.45,      0.30,       0.02,
  "Recovery 1",            1.00,      0.10,       0.02,            0.80,      0.10,       0.40,
  "Treatment 2",           0.45,      0.47,       0.30,            0.45,      0.30,       0.02,
  "Recovery 2",            1.00,      0.10,       0.02,            0.80,      0.10,       0.40
)

# -----------------------------------------------------------------------------
# 3. Initialize tissue
# -----------------------------------------------------------------------------
sim <- TissueSimulation(
  epigenetic_states = c("P1", "P2"),
  width = 2000,
  height = 2000,
  save_directory = TRUE,
  name = "epigenetic_memory_double_treatment"
)
sim$death_activation_level <- 50
sim$add_mutant("G1", rates_baseline)
sim$place_cell("G1[P1]", 1000, 1000)

# -----------------------------------------------------------------------------
# 4. Baseline
# -----------------------------------------------------------------------------
count_initial <- get_state_counts(sim, "Baseline")
count_baseline <- run_and_record(sim, TIMELINE$baseline_end, "Baseline", step = record_every)
p_pre <- plot_state(sim, color_map = STATE_COLORS) + labs(title = "Pre-treatment")
sample_pre <- sample_tissue_plot(sim, "S_PRE_T1", STATE_COLORS)

# -----------------------------------------------------------------------------
# 5. Treatment 1
# -----------------------------------------------------------------------------
sim$set_rates(list(G1 = rates_T1))
count_T1_start <- get_state_counts(sim, "Treatment 1", 1)
count_T1 <- run_and_record(sim, TIMELINE$T1_end, "Treatment 1", 1, record_every)
p_T1 <- plot_state(sim, color_map = STATE_COLORS) + labs(title = "Post-Treatment 1")
sample_post_T1 <- sample_tissue_plot(sim, "S_POST_T1", STATE_COLORS)

# -----------------------------------------------------------------------------
# 6. Recovery 1 / relapse
# -----------------------------------------------------------------------------
sim$set_rates(list(G1 = rates_recovery))
count_recovery1 <- run_and_record(sim, TIMELINE$recovery1_end, "Recovery 1", step = record_every)
p_pre_T2 <- plot_state(sim, color_map = STATE_COLORS) + labs(title = "Pre-Treatment 2 / relapse")
sample_pre_T2 <- sample_tissue_plot(sim, "S_PRE_T2", STATE_COLORS)

# -----------------------------------------------------------------------------
# 7. Treatment 2: same fitness rates, increased P1 -> P2 switching
# -----------------------------------------------------------------------------
sim$set_rates(list(G1 = rates_T2))
count_T2_start <- get_state_counts(sim, "Treatment 2", 2)
count_T2 <- run_and_record(sim, TIMELINE$T2_end, "Treatment 2", 2, record_every)
p_T2 <- plot_state(sim, color_map = STATE_COLORS) + labs(title = "Post-Treatment 2")
sample_post_T2 <- sample_tissue_plot(sim, "S_POST_T2", STATE_COLORS)

# -----------------------------------------------------------------------------
# 8. Recovery 2
# -----------------------------------------------------------------------------
sim$set_rates(list(G1 = rates_recovery))
count_recovery2 <- run_and_record(sim, TIMELINE$recovery2_end, "Recovery 2", step = record_every)
p_final <- plot_state(sim, color_map = STATE_COLORS) + labs(title = "Second relapse")
sample_final <- sample_tissue_plot(sim, "S_RELAPSE_T2", STATE_COLORS)
sample_forest <- sim$get_sample_forest()

# -----------------------------------------------------------------------------
# 9. Assemble results
# -----------------------------------------------------------------------------
count_time <- bind_rows(
  count_initial, count_baseline,
  count_T1_start, count_T1,
  count_recovery1,
  count_T2_start, count_T2,
  count_recovery2
) |>
  mutate(
    time = as.numeric(time),
    epistate = factor(epistate, levels = c("P1", "P2"))
  ) |>
  arrange(time, epistate)

treatment_dynamics <- prepare_treatment_dynamics(count_time, TIMELINE)
rates_history <- prepare_rate_history(sim)

# -----------------------------------------------------------------------------
# 10. Plots
# -----------------------------------------------------------------------------
max_Y <- max(sample_forest$get_nodes()$birth_time, na.rm = TRUE)
to_y  <- function(t) max_Y - t   # time -> plot coordinate

plot_sample_forest <-plot_forest(sample_forest, color_map = STATE_COLORS,highlight_sample = T) +
  # Treatment 1
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = to_y(TIMELINE$T1_end),
           ymax = to_y(TIMELINE$T1_start),
           alpha = 0.2) +
  # Treatment 2
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = to_y(TIMELINE$T2_end),
           ymax = to_y(TIMELINE$T2_start),
           alpha = 0.2)
plot_sample_forest <- annotate_forest(plot_sample_forest,forest = sample_forest,samples = T,MRCAs = T)

p_population <- plot_total_population(count_time, TIMELINE)
p_absolute <- plot_absolute_states(count_time, TIMELINE, EPISTATE_COLORS)
p_composition <- plot_composition(count_time, TIMELINE, EPISTATE_COLORS)
p_P2_memory <- plot_p2_memory(treatment_dynamics, treatment_duration)
p_memory_delta <- plot_p2_delta(treatment_dynamics, treatment_duration)
p_normalized_killing <- plot_normalized_burden(treatment_dynamics, treatment_duration)
p_rates_growth <- plot_growth_rates(rates_history, TIMELINE, EPISTATE_COLORS)
p_rates_switch <- plot_switch_rates(rates_history, TIMELINE)

p_spatial <- p_pre | p_T1 | p_pre_T2 | p_T2 | p_final
p_samples <- sample_pre | sample_post_T1 | sample_pre_T2 | sample_post_T2 | sample_final

p_overview <- (p_population / p_absolute / p_composition) +
  plot_annotation(title = "Repeated-treatment simulation overview",tag_levels = "a")

p_treatment_comparison <- (p_normalized_killing / p_P2_memory / p_memory_delta) +
  plot_annotation(title = "Treatment 1 vs Treatment 2",tag_levels = "a")

p_rate_summary <- (p_rates_growth / p_rates_switch) +
  plot_annotation(title = "Model parameters over time",tag_levels = "a")

# -----------------------------------------------------------------------------
# 11. Save reproducible outputs
# -----------------------------------------------------------------------------
dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create(file.path("results", "figures"), showWarnings = FALSE, recursive = TRUE)

readr::write_csv(count_time, file.path("results", "count_time.csv"))
readr::write_csv(parameter_summary, file.path("results", "parameter_summary.csv"))
readr::write_csv(rates_history, file.path("results", "rates_history.csv"))
saveRDS(sim, file.path("results", "simulation.rds"))
saveRDS(count_time, file.path("results", "count_time.rds"))

save_plot(p_overview, file.path("results", "figures", "01_overview.pdf"), 10, 12)
save_plot(p_treatment_comparison, file.path("results", "figures", "02_T1_vs_T2.pdf"), 9, 11)
save_plot(p_rate_summary, file.path("results", "figures", "03_rate_history.pdf"), 10, 8)
save_plot(p_spatial, file.path("results", "figures", "04_spatial_evolution.pdf"), 16, 4.5)
save_plot(p_samples, file.path("results", "figures", "05_sampled_tissues.pdf"), 16, 4.5)
save_plot(plot_sample_forest, file.path("results", "figures", "06_sample_forest.pdf"),5, 6)

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      # Save all objects required by the report in one file.
save(
  count_time, treatment_dynamics, rates_history, parameter_summary,
  p_population, p_absolute, p_composition,
  p_P2_memory, p_memory_delta, p_normalized_killing,
  p_rates_growth, p_rates_switch, p_spatial, p_samples,
  p_overview, p_treatment_comparison, p_rate_summary,
  plot_sample_forest,
  file = file.path("results", "report_objects.RData")
)
sample_forest$save(file.path(prj_dir,"results","sample_forest.sff"))
message("Simulation complete. Results written to ./results/")
message("Render final_report.Rmd after this script finishes.")




##################### MUTATION PART ##############
setwd("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project/ProCESS-scATAC/scripts/simulate_fragments/peak_based/0_process_simulations/process_references_v1.3.5/")

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


## Drivers for the tumors
m_engine$add_mutant(mutant_name = "G1",
                    passenger_rates = list("P1" = c(SNV = mu_SNV, indel = mu_INDELs,CNA=mu_CNA),
                                           "P2" = c(SNV = mu_SNV, indel = mu_INDELs,CNA=mu_CNA)),
                    drivers = list("APC Q1294Gfs*6","TP53 R175H",CNA_Clone0_1,CNA_Clone0_2,CNA_Clone0_3))
m_engine$change_rates_from(TIMELINE$T1_start, "G1", list("P2" = c(SNV = 5e-9, indel = 2e-9)))
m_engine$change_rates_from(TIMELINE$T1_end, "G1", list("P2" = c(SNV = mu_SNV, indel = mu_INDELs)))
m_engine$change_rates_from(TIMELINE$T2_start, "G1", list("P2" = c(SNV = 5e-9, indel = 4e-9)))
m_engine$change_rates_from(TIMELINE$T2_end, "G1", list("P2" = c(SNV = mu_SNV, indel = mu_INDELs)))
# ============================================================
# PRE-TREATMENT: baseline colorectal cancer
# ============================================================

m_engine$add_exposure(
  time = 0,
  coefficients = c(
    
    # SBS signatures — sum = 1
    SBS1  = 0.30,   # clock-like
    SBS5  = 0.35,   # clock-like
    SBS18 = 0.20,   # oxidative damage / ROS
    SBS88 = 0.15,   # colibactin
    
    # ID signatures — sum = 1
    ID1   = 0.40,   # replication slippage
    ID2   = 0.35,   # replication slippage
    ID18  = 0.25    # colibactin
  )
)


# ============================================================
# DURING FOLFOX
# 5-FU + oxaliplatin
# ============================================================

m_engine$add_exposure(
  time = TIMELINE$T1_start,
  coefficients = c(
    
    # SBS signatures — sum = 1
    SBS1   = 0.10,
    SBS5   = 0.15,
    SBS18  = 0.10,
    SBS88  = 0.05,
    
    SBS17b = 0.35,  # 5-FU-associated
    SBS35  = 0.25,  # platinum-associated
    
    # ID signatures — sum = 1
    ID1    = 0.40,
    ID2    = 0.35,
    ID18   = 0.25
  )
)


# ============================================================
# POST-TREATMENT T1
# chemotherapy-associated mutational processes switched off
# ============================================================

m_engine$add_exposure(
  time = TIMELINE$T1_end,
  coefficients = c(
    
    # SBS signatures — sum = 1
    SBS1   = 0.30,
    SBS5   = 0.35,
    SBS18  = 0.20,
    SBS88  = 0.15,
    
    # Treatment processes inactive
    SBS17b = 0.00,
    SBS35  = 0.00,
    
    # ID signatures — sum = 1
    ID1    = 0.40,
    ID2    = 0.35,
    ID18   = 0.25
  )
)

# ============================================================
# DURING FOLFOX
# 5-FU + oxaliplatin
# ============================================================

m_engine$add_exposure(
  time = TIMELINE$T2_start,
  coefficients = c(
    
    # SBS signatures — sum = 1
    SBS1   = 0.1,
    SBS5   = 0.1,
    SBS18  = 0.05,
    SBS88  = 0.05,
    
    SBS17b = 0.40,  # 5-FU-associated
    SBS35  = 0.30,  # platinum-associated
    
    # ID signatures — sum = 1
    ID1    = 0.40,
    ID2    = 0.35,
    ID18   = 0.25
  )
)


# ============================================================
# POST-TREATMENT T1
# chemotherapy-associated mutational processes switched off
# ============================================================

m_engine$add_exposure(
  time = TIMELINE$T2_end,
  coefficients = c(
    
    # SBS signatures — sum = 1
    SBS1   = 0.30,
    SBS5   = 0.35,
    SBS18  = 0.20,
    SBS88  = 0.15,
    
    # Treatment processes inactive
    SBS17b = 0.00,
    SBS35  = 0.00,
    
    # ID signatures — sum = 1
    ID1    = 0.40,
    ID2    = 0.35,
    ID18   = 0.25
  )
)

if (!exists("sample_forest")){
  sample_forest <- load_sample_forest(file.path(prj_dir,"results","sample_forest.sff"))
}
phylo_forest <- m_engine$place_mutations(sample_forest, num_of_preneoplatic_SNVs=800, num_of_preneoplatic_indels=200)

setwd(prj_dir)
phylo_forest$save(file.path(prj_dir,"results","phylo_forest.sff"))

############# PLOTS #############
plot_exposure <- plot_exposure_timeline(
  phylo_forest,
  linewidth = 1.5,
  emphasize_switches = TRUE,
  mutation_type = "SNV",
  pal_name = "Set3"
)+
  annotate(
    "rect",
    xmin = TIMELINE$T1_start,
    xmax = TIMELINE$T1_end,
    ymin = -Inf,
    ymax = Inf,
    alpha = 0.07
  )+
  annotate(
    "rect",
    xmin = TIMELINE$T2_start,
    xmax = TIMELINE$T2_end,
    ymin = -Inf,
    ymax = Inf,
    alpha = 0.07
  )

species_info <- m_engine$get_species_info()
mut_rate_df <- species_info %>%
  select(mutant, epistate, time, SNV_rate, indel_rate) %>%
  pivot_longer(
    cols = c(SNV_rate, indel_rate),
    names_to = "mutation_type",
    values_to = "mutation_rate"
  ) %>%
  dplyr::mutate(mutation_rate_log=log10(mutation_rate)) %>% 
  mutate(
    mutation_type = recode(
      mutation_type,
      "SNV_rate"   = "SNV",
      "indel_rate" = "Indel"
    )
  ) %>% 
  group_by(epistate, mutation_type) %>%
  group_modify(~ bind_rows(
    .x,
    .x %>%
      slice_max(time, n = 1) %>%
      mutate(time = 200)
  )) %>%
  ungroup()
last_mut_rate_df <- mut_rate_df %>%
  group_by(mutation_type) %>%
  slice_max(time, n = 1, with_ties = FALSE) %>%
  mutate(time = sim$get_clock())

mut_rate_df <- bind_rows(
  mut_rate_df,
  last_mut_rate_df
)
p_mut_rate_hist <-ggplot(
  mut_rate_df,
  aes(
    x = time,
    y = mutation_rate,
    color = epistate,
    group = epistate
  )
) +
  
  geom_step(
    direction = "hv",
    linewidth = 1.2
  ) +
  
  facet_wrap(
    ~ mutation_type,
    ncol = 1,
    scales = "free_y"
  ) +
  
  scale_color_manual(
    values = c(
      "P1" = "goldenrod",
      "P2" = "orchid2"
    )
  ) +
  annotate(
    "rect",
    xmin = TIMELINE$T1_start,
    xmax = TIMELINE$T1_end,
    ymin = -Inf,
    ymax = Inf,
    alpha = 0.07
  )+
  annotate(
    "rect",
    xmin = TIMELINE$T2_start,
    xmax = TIMELINE$T2_end,
    ymin = -Inf,
    ymax = Inf,
    alpha = 0.07
  )+
  labs(
    x = "Simulation time",
    y = "Mutation rate",
    color = "Phenotype"
  ) +
  theme_minimal()+
  
  theme(
    legend.position = "right",
  )
p_genomics <- (p_mut_rate_hist / plot_exposure) +
  plot_annotation(title = "Genomic-paramters evolution",tag_levels = "a")

save_plot(p_genomics, file.path(prj_dir,"results", "figures", "07_genomic_evo.pdf"), 10, 8)

