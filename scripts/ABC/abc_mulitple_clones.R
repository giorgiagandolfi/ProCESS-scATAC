rm(list=ls())
library(ProCESS)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(EasyABC)

# create a tissue simulation with two epigenetic states
sim <- TissueSimulation(width = 1000, height = 1000,save_snapshots = F,name = "P05")

rates_genetic_clone_0 <- list(duplication = 0.6, death = 0.1)


# add a mutant "A" and set its species rates
rates_genetic_clone_A <- list(
  duplication = 1.5, death = 0.1
)
sim$add_mutant("0", rates_genetic_clone_0)

sim$place_cell("0", 500, 500)

# let the simulation evolve until the species "A[E2]" has less than 10 cells
sim$run_up_to_time(50)


sim$add_mutant("A",rates_genetic_clone_A)
starting_cell=sim$choose_border_cell_in("0")
plot_tissue(sim)+
  # geom_vline(xintercept = starting_cell$position_x)+
  # geom_hline(yintercept = starting_cell$position_y)+
  geom_point(data = starting_cell,aes(x=position_x,y=position_y),inherit.aes = FALSE)

sim$mutate_progeny(starting_cell,"A")
sim$run_up_to_time(120)
plot_tissue(sim)


b_A=rates_genetic_clone_A$duplication
d_A=rates_genetic_clone_A$death
tot_A=sim$get_count_history() %>% 
  filter(mutant=="A") %>% 
  tail(1) %>% pull(count)
t_A=120.00002-(log(tot_A)/(b_A-d_A))
t_0=50.00014-(log(11943)/(rates_genetic_clone_0$duplication-rates_genetic_clone_0$death))

# Sample Clone A
n_w <- n_h <- 20
ncells <- 0.8 * n_w * n_h
bbox0 <- sim$search_sample(c("0" = ncells), n_w, n_h)
bboxA <- sim$search_sample(c("A" = ncells), n_w, n_h)
sim$sample_cells("S_0", bbox0$lower_corner, bbox0$upper_corner)
sim$sample_cells("S_A", bboxA$lower_corner, bboxA$upper_corner)
plot_tissue(sim)

sample_forest <- sim$get_sample_forest()
observed=sample_forest$get_nodes() %>% filter(!is.na(sample)) %>% count(mutant)



proccess_model <- function(x) {
  print(x)
  # x contains:
  # x[1] = duplication_0
  # x[2] = death_0
  # x[3] = duplication_A
  # x[4] = death_A
  
  duplication_0 <- x[2]
  death_0       <- x[3]
  
  duplication_A <- x[4]
  death_A       <- x[5]
  
  print(duplication_0)
  print(death_0)
  print(duplication_A)
  print(death_A)
  # ------------------------------------
  # Create simulation
  # ------------------------------------
  
  sim <- TissueSimulation(
    width = 1000,
    height = 1000,
    save_snapshots = FALSE,
    name = paste0("ABC_", Sys.getpid())
  )
  
  
  # ------------------------------------
  # Define rates
  # ------------------------------------
  
  rates_0 <- list(
    duplication = duplication_0,
    death = death_0
  )
  
  rates_A <- list(
    duplication = duplication_A,
    death = death_A
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
  
  sim$run_up_to_time(10)
  
  
  # ------------------------------------
  # Introduce clone A
  # ------------------------------------
  
  sim$add_mutant(
    "A",
    rates_A
  )
  
  starting_cell <- sim$choose_border_cell_in("0")
  
  sim$mutate_progeny(
    starting_cell,
    "A"
  )
  
  
  # ------------------------------------
  # Continue evolution
  # ------------------------------------
  
  sim$run_up_to_time(30)
  
  
  # ------------------------------------
  # Sample
  # ------------------------------------
  
  # Sample Clone A
  n_w <- n_h <- 5
  ncells <- 0.8 * n_w * n_h
  bbox0 <- sim$search_sample(c("0" = ncells), n_w, n_h)
  bboxA <- sim$search_sample(c("A" = ncells), n_w, n_h)
  sim$sample_cells("S_0", bbox0$lower_corner, bbox0$upper_corner)
  sim$sample_cells("S_A", bboxA$lower_corner, bboxA$upper_corner)
  
  
  sample_forest <- sim$get_sample_forest()
  
  
  
  # ------------------------------------
  # Get sampled cells
  # ------------------------------------
  
  sample_forest <- sim$get_sample_forest()
  
  counts <- sample_forest$get_nodes() %>%
    filter(!is.na(sample)) %>%
    count(mutant)
  
  
  # ------------------------------------
  # Make sure both clones exist
  # ------------------------------------
  
  n_0 <- counts %>%
    filter(mutant == "0") %>%
    pull(n)
  
  n_A <- counts %>%
    filter(mutant == "A") %>%
    pull(n)
  
  
  if (length(n_0) == 0)
    n_0 <- 0
  
  if (length(n_A) == 0)
    n_A <- 0
  
  
  # ------------------------------------
  # Return ABC summary statistics
  # ------------------------------------
  
  c(
    n_0 = n_0,
    n_A = n_A
  )
}
invalid_summary <- function() {
  c(
    n_0 = 1e6,
    n_A = 1e6
  )
}
proccess_model_trycatch <- function(x) {
  
  # x[1] = seed
  # x[2] = duplication_0
  # x[3] = death_0
  # x[4] = duplication_A
  # x[5] = death_A
  
  set.seed(x[1])
  
  duplication_0 <- x[2]
  death_0       <- x[3]
  duplication_A <- x[4]
  death_A       <- x[5]
  
  # ------------------------------------
  # Create simulation
  # ------------------------------------
  
  sim <- TissueSimulation(
    width = 1000,
    height = 1000,
    save_snapshots = FALSE,
    name = paste0("ABC_", Sys.getpid())
  )
  
  rates_0 <- list(
    duplication = duplication_0,
    death = death_0
  )
  
  rates_A <- list(
    duplication = duplication_A,
    death = death_A
  )
  
  # ------------------------------------
  # Clone 0
  # ------------------------------------
  
  sim$add_mutant("0", rates_0)
  
  sim$place_cell("0", 500, 500)
  
  sim$run_up_to_time(50)
  
  
  # ------------------------------------
  # Find border cell
  # ------------------------------------
  
  starting_cell <- tryCatch(
    
    sim$choose_border_cell_in("0"),
    
    error = function(e) {
      NULL
    }
  )
  
  # No suitable cell -> reject proposal
  
  if (is.null(starting_cell)) {
    return(invalid_summary())
  }
  
  
  # ------------------------------------
  # Introduce clone A
  # ------------------------------------
  
  sim$add_mutant("A", rates_A)
  
  mutation_result <- tryCatch(
    
    sim$mutate_progeny(
      starting_cell,
      "A"
    ),
    
    error = function(e) {
      NULL
    }
  )
  
  if (is.null(mutation_result)) {
    return(invalid_summary())
  }
  
  
  # ------------------------------------
  # Continue evolution
  # ------------------------------------
  
  sim$run_up_to_time(120)
  
  
  # ------------------------------------
  # Find sampling boxes
  # ------------------------------------
  
  n_w <- 20
  n_h <- 20
  ncells <- 0.8 * n_w * n_h
  
  bbox0 <- tryCatch(
    
    sim$search_sample(
      c("0" = ncells),
      n_w,
      n_h
    ),
    
    error = function(e) {
      NULL
    }
  )
  
  if (is.null(bbox0)) {
    return(invalid_summary())
  }
  
  
  bboxA <- tryCatch(
    
    sim$search_sample(
      c("A" = ncells),
      n_w,
      n_h
    ),
    
    error = function(e) {
      NULL
    }
  )
  
  if (is.null(bboxA)) {
    return(invalid_summary())
  }
  
  
  # ------------------------------------
  # Sample
  # ------------------------------------
  
  sim$sample_cells(
    "S_0",
    bbox0$lower_corner,
    bbox0$upper_corner
  )
  
  sim$sample_cells(
    "S_A",
    bboxA$lower_corner,
    bboxA$upper_corner
  )
  
  
  # ------------------------------------
  # Get counts
  # ------------------------------------
  
  sample_forest <- sim$get_sample_forest()
  
  counts <- sample_forest$get_nodes() %>%
    filter(!is.na(sample)) %>%
    count(mutant)
  
  n_0 <- counts %>%
    filter(mutant == "0") %>%
    pull(n)
  
  n_A <- counts %>%
    filter(mutant == "A") %>%
    pull(n)
  
  
  if (length(n_0) == 0) n_0 <- 0
  if (length(n_A) == 0) n_A <- 0
  
  
  # ------------------------------------
  # Return summary
  # ------------------------------------
  
  c(
    n_0 = n_0,
    n_A = n_A
  )
}
prior <- list(
  c("unif", 0.5, 0.8),   # duplication_0
  c("unif", 0.08, 0.2),  # death_0
  
  c("unif", 1, 1.7),   # duplication_A
  c("unif", 0.08, 0.2)   # death_A
)

summary_stat_target <- c(
  n_0 = observed$n[1],
  n_A = observed$n[2]
)

ABC_result <- ABC_sequential(
  method = "Beaumont",
  model = proccess_model_trycatch,
  prior = prior,
  nb_simul = 50,
  summary_stat_target = summary_stat_target,
  tolerance_tab = c(
    50,10
  ),
  use_seed = TRUE,
  n_cluster = 1,
  verbose = TRUE
)


posterior <- as.data.frame(ABC_result$param)

colnames(posterior) <- c(
  "duplication_0",
  "death_0",
  "duplication_A",
  "death_A"
)

posterior_long <- posterior %>%
  pivot_longer(
    cols = everything(),
    names_to = "parameter",
    values_to = "value"
  )
ggplot(
  posterior_long,
  aes(x = value)
) +
  geom_histogram(bins = 30) +
  facet_wrap(
    ~ parameter,
    scales = "free"
  ) +
  theme_bw()



prior_df <- tibble(
  parameter = c(
    "duplication_0",
    "death_0",
    "duplication_A",
    "death_A"
  ),
  true_param = c(rates_genetic_clone_0$duplication,rates_genetic_clone_0$death,rates_genetic_clone_A$duplication,rates_genetic_clone_A$death)
)
ggplot(posterior_long, aes(x = value)) +
  geom_histogram(bins = 15) +
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
