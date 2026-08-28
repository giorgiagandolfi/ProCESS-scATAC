rm(list=ls())
library(ProCESS)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(EasyABC)
setwd("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project/ProCESS-scATAC/scripts/ABC/")
# create a tissue simulation with two epigenetic states
# sim <- TissueSimulation(width = 1000, height = 1000,save_snapshots = F,name = "P05")
sim <- TissueSimulation(width = 1000, height = 1000,name = "P05")

rates_genetic_clone_0 <- list(duplication = 0.6, death = 0.1)
sim$add_mutant("0", rates_genetic_clone_0)

sim$place_cell("0", 500, 500)

# let the simulation evolve until the species "A[E2]" has less than 10 cells
sim$run_up_to_time(50)

tot_counts_observed=sim$get_count_history()
plot_tissue(sim)


proccess_model_small <- function(x) {
  library(ProCESS)
  library(dplyr)
  print(x)
  duplication_0 <- x[2]
  death_0       <- x[3]

  print(duplication_0)
  print(death_0)

  sim <- TissueSimulation(
    width = 1000,
    height = 1000,
    name = paste0("ABC_", Sys.getpid())
  )
  
  
  # ------------------------------------
  # Define rates
  # ------------------------------------
  
  rates_0 <- list(
    duplication = duplication_0,
    death = death_0
  )
  

  
  # ------------------------------------
  # Clone 0
  # ------------------------------------
  
  sim$add_mutant(
    "0",
    rates_0
  )
  
  sim$place_cell(
    "0",
    500,
    500
  )
  
  
  # ------------------------------------
  # Grow clone 0
  # ------------------------------------
  
  sim$run_up_to_time(50)
  

  
  
  
  # ------------------------------------
  # Make sure both clones exist
  # ------------------------------------
  tot_counts1=sim$get_count_history()
  n_0 <- tot_counts1 %>%
    filter(mutant == "0") %>%
    pull(count)
  

  
  
  if (length(n_0) == 0)
    n_0 <- 0

  # ------------------------------------
  # Return ABC summary statistics
  # ------------------------------------
  
  c(
    n_0 = n_0
  )
}
prior <- list(
  c("unif", 0.3, 2),   # duplication_0
  c("unif", 0.01, 0.3)  # death_0
)

summary_stat_target <- c(
  n_0 = tot_counts_observed$count[1]
)

# ABC_result <- ABC_sequential(
#   method = "Beaumont",
#   model = proccess_model_small,
#   prior = prior,
#   nb_simul = 1e3,
#   summary_stat_target = summary_stat_target,
#   tolerance_tab = c(2, 1, 0.5, 0.25),
#   use_seed = TRUE,
#   n_cluster = 1,
#   verbose = TRUE
# )

ABC_rej<-ABC_rejection(model=proccess_model_small, prior=prior, nb_simul=1e4,
                       summary_stat_target=summary_stat_target, tol=.3,use_seed = T)
posterior <- as.data.frame(ABC_rej$param)

colnames(posterior) <- c(
  "duplication_0",
  "death_0"
)

posterior_long <- posterior %>%
  pivot_longer(
    cols = everything(),
    names_to = "parameter",
    values_to = "value"
  )



prior_df <- tibble(
  parameter = c(
    "duplication_0",
    "death_0"
  ),
  true_param = c(rates_genetic_clone_0$duplication,rates_genetic_clone_0$death)
)

p=ggplot(posterior_long, aes(x = value)) +
  geom_histogram() +
  geom_vline(
    data = prior_df,
    aes(xintercept = true_param),
    linetype = "dashed"
  ) +
  facet_wrap(
    ~ parameter,
    scales = "free"
  ) +
  theme_bw()
ggsave(filename = "abc_single_clone_rejection.pdf",plot = p)