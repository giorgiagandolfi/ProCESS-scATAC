rm(list=ls())
library(ProCESS)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggalluvial)
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project/ProCESS-scATAC/scripts/simulate_fragments/peak_based/utils.R")
set.seed(0)

sim <- TissueSimulation(epigenetic_states = c("E1", "E2","E3"),
                        width = 1000, height = 1000,save_snapshots = F,name = "case2")
sim$death_activation_level <- 50
rates_genetic_clone_0 <- list(
  E1 = list(duplication = 0.8, death = 0.1, E2 = 0.8, E3 = 0.7),
  E2 = list(duplication = 0.8, death = 0.1, E1 = 0.7, E3 = 0.8),
  E3 = list(duplication = 0.8, death = 0.1, E1 = 0.8, E2 = 0.7)
)

# clone A: low switching, single major epigenetic state = E1
rates_genetic_clone_A <- list(
  E1 = list(duplication = 2, death = 0.1, E2 = 0.02, E3 = 0.02),
  E2 = list(duplication = 2, death = 0.1, E1 = 0.5,  E3 = 0.05),
  E3 = list(duplication = 2, death = 0.1, E1 = 0.5,  E2 = 0.05)
)

# clone B: low switching, single major epigenetic state = E3 (different from A)
rates_genetic_clone_B <- list(
  E1 = list(duplication = 2.5, death = 0.1, E2 = 0.05, E3 = 0.5),
  E2 = list(duplication = 2.5, death = 0.1, E1 = 0.05, E3 = 0.5),
  E3 = list(duplication = 2.5, death = 0.1, E1 = 0.02, E2 = 0.02)
)

q_G1=get_epigenetic_Q(epigenetic_rates = rates_genetic_clone_0,epistate_name = "G1")
q_G2=get_epigenetic_Q(epigenetic_rates = rates_genetic_clone_A,epistate_name = "G2")
q_G3=get_epigenetic_Q(epigenetic_rates = rates_genetic_clone_B,epistate_name = "G3")
q_G1$steady_states
q_G2$steady_states
q_G3$steady_states


sim$add_mutant("G1", rates_genetic_clone_0)
sim$place_cell("G1[E2]", 500, 500)

# let the simulation evolve until the species "A[E2]" has less than 10 cells
sim$run_up_to_time(10)
sim$get_counts()

sim$add_mutant("G2",rates_genetic_clone_A)
starting_cell=sim$choose_border_cell_in("G1[E1]")
plot_tissue(sim)+
  geom_point(data = starting_cell,aes(x=position_x,y=position_y),inherit.aes = FALSE)+
  facet_wrap(~epistate)

sim$mutate_progeny(starting_cell,"G2")
sim$set_rates(list("G1" = list(duplication = 0.4,death=0.2)))
sim$run_up_to_time(40)
plot_tissue(sim)+
  facet_wrap(~epistate)
sim$get_counts()

sim$add_mutant("G3",rates_genetic_clone_B)
starting_cell=sim$choose_border_cell_in("G2[E3]")
plot_tissue(sim)+
  geom_point(data = starting_cell,aes(x=position_x,y=position_y),inherit.aes = FALSE)+
  facet_wrap(~epistate)

sim$mutate_progeny(starting_cell,"G3")
sim$set_rates(list("G2" = list(duplication = 0.4, death = 0.2)))
sim$set_rates(list("G1" = list(duplication = 0.3, death = 0.3)))
sim$run_up_to_time(60)
plot_tissue(sim)+
  facet_grid(mutant~epistate)



# # Sample Clone 0
# n_w <- n_h <- 15
# ncells <- 0.8 * n_w * n_h
# bbox <- sim$search_sample(c("0" = ncells), n_w, n_h)
# sim$sample_cells("S1", bbox$lower_corner, bbox$upper_corner)
# plot_tissue(sim,at_sample = "S1")
# 
# # Sample Clone B
# n_w <- n_h <- 15
# ncells <- 0.8 * n_w * n_h
# bbox <- sim$search_sample(c("A" = ncells), n_w, n_h)
# sim$sample_cells("S2", bbox$lower_corner, bbox$upper_corner)
# plot_tissue(sim,at_sample = "S2")
# 
# # Sample Clone B
# n_w <- n_h <- 15
# ncells <- 0.8 * n_w * n_h
# bbox <- sim$search_sample(c("B" = ncells), n_w, n_h)
# sim$sample_cells("S3", bbox$lower_corner, bbox$upper_corner)
# plot_tissue(sim,at_sample = "S3")

# Sample all clones toghter
bbox_width=40
bbox1_p <- c(430,405)
bbox1_q <- bbox1_p + bbox_width
# view the boxes
plot_tissue(sim) +
  geom_rect(
    xmin = bbox1_p[1],
    xmax = bbox1_q[1],
    ymin = bbox1_p[2],
    ymax = bbox1_q[2],
    fill = NA,
    color = "black"
  )
sim$sample_cells("S3", bbox1_p, bbox1_q)
plot_tissue(sim,at_sample = "S3")
sample_forest <- sim$get_sample_forest()

