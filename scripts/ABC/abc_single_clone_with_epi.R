rm(list=ls())
library(ProCESS)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(EasyABC)
setwd("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project/ProCESS-scATAC/scripts/ABC/")
# create a tissue simulation with two epigenetic states
sim <- TissueSimulation(epigenetic_states = c("E1", "E2"),width = 1000, height = 1000,name = "P05")
# sim <- TissueSimulation(width = 1000, height = 1000,save_directory=F,name = "P05")


rates_genetic_clone_0 <- list(
  E1 = list(duplication = 1.2, death = 0.1, E2 = 0.8),
  E2 = list(duplication = 1.2, death = 0.1, E1 = 0.3)
  # E3 = list(duplication = 0.8, death = 0.1, E1 = 0.01, E2 = 0.2, E4 = 0.02),
  # E4 = list(duplication = 0.8, death = 0.1, E1 = 0.02, E2 = 0.2, E3 = 0.5)
)
sim$add_mutant("0", rates_genetic_clone_0)

sim$place_cell("0[E1]", 500, 500)

# let the simulation evolve until the species "A[E2]" has less than 10 cells
sim$run_up_to_time(50)
plot_tissue(sim)+facet_wrap(~epistate)
tot_counts_observed=sim$get_count_history()
tot_counts_observed
proccess_model_small_epi <- function(x) {
  library(ProCESS)
  library(dplyr)

  print(x)
  switch_e1_e2 <- x[2]
  switch_e2_e1       <- x[3]

  
  sim <- TissueSimulation(
    epigenetic_states = c("E1", "E2"),
    width = 1000,
    height = 1000,
    name = paste0("ABC_", Sys.getpid())
  )
  
  
  # ------------------------------------
  # Define rates
  # ------------------------------------
  
  rates_genetic_clone_0 <- list(
    E1 = list(duplication = 1.2, death = 0.1, E2 = switch_e1_e2),
    E2 = list(duplication = 1.2, death = 0.1, E1 = switch_e2_e1)
  )
  
  
  
  # ------------------------------------
  # Clone 0
  # ------------------------------------
  
  sim$add_mutant(
    "0",
    rates_genetic_clone_0
  )
  
  sim$place_cell("0[E1]", 500, 500)
  
  sim$run_up_to_time(50)
  
  
  
  
  
  # ------------------------------------
  # Make sure both clones exist
  # ------------------------------------
  tot_counts_epi=sim$get_count_history()
  n_e1 <- tot_counts_epi %>%
    filter(epistate == "E1") %>%
    pull(count)
  n_e2 <- tot_counts_epi %>%
    filter(epistate == "E2") %>%
    pull(count)
  
  
  
  
  if (length(n_e1) == 0)
    n_e1 <- 0
  if (length(n_e2) == 0)
    n_e2 <- 0
  
  # ------------------------------------
  # Return ABC summary statistics
  # ------------------------------------
  print(n_e1)
  print(n_e2)
  c(
    n_e1 = n_e1,
    n_e2 = n_e2
  )
}
prior <- list(
  c("unif", 0.6, 1),   # duplication_0
  c("unif", 0.2, 0.4)  # death_0
)

summary_stat_target <- c(
  n_e1 = tot_counts_observed$count[1],
  n_e2 = tot_counts_observed$count[2]
)
ABC_result <- ABC_sequential(
  method = "Beaumont",
  model = proccess_model_small_epi,
  prior = prior,
  nb_simul = 100,
  summary_stat_target = summary_stat_target,
  tolerance_tab = c(2, 1, 0.5, 0.25),
  use_seed = TRUE,
  n_cluster = 1,
  verbose = TRUE
)
posterior <- as.data.frame(ABC_result$param)

colnames(posterior) <- c(
  "e1_e2",
  "e2_e1"
)

posterior_long <- posterior %>%
  pivot_longer(
    cols = everything(),
    names_to = "parameter",
    values_to = "value"
  )



prior_df <- tibble(
  parameter = c(
    "e1_e2",
    "e2_e1"
  ),
  true_param = c(rates_genetic_clone_0$E1$E2,rates_genetic_clone_0$E2$E1)
)

p=ggplot(posterior_long, aes(x = value)) +
  geom_histogram() +
  geom_vline(
    data = prior_df,
    aes(xintercept = true_param),
    linetype = "dashed",colour = "red"
  ) +
  facet_wrap(
    ~ parameter,
    scales = "free"
  ) +
  theme_bw()
ggsave(filename = "abc_single_clone_two_epistate.pdf",plot = p)