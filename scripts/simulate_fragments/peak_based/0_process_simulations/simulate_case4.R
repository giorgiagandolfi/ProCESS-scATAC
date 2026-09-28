rm(list=ls())
library(ProCESS)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggalluvial)
library(patchwork)
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project/ProCESS-scATAC/scripts/simulate_fragments/peak_based/utils.R")
set.seed(0)

sim <- TissueSimulation(epigenetic_states = c("P1", "P2"),
                        width = 2000, height = 2000,save_directory = T,name = "case4")
sim$death_activation_level <- 50
rates_genetic_clone_G1 <- list(
  P1 = list(duplication = 1.2, death = 0.1, P2 = 0.01),
  P2 = list(duplication = 1.2, death = 0.1, P1 = 0.9)
)

sim$add_mutant("G1", rates_genetic_clone_G1)
sim$place_cell("G1[P1]", 1000, 1000)
sim$run_up_to_time(50)

p1=plot_state(sim,color_map = c("G1[P1]"="goldenrod","G1[P2]"="orchid2"))
count_time_p1 <- sim$get_cells() %>% 
  dplyr::group_by(epistate) %>% 
  dplyr::summarise(n=n()) %>% 
  dplyr::mutate(prop = n / sum(n)) %>% 
  dplyr::mutate(time=sim$get_clock())


sim$set_rates(list("G1" = list(duplication = 0.8, death = 0.3)))
sim$run_up_to_time(80)
sim$get_counts()
p2=plot_state(sim,color_map = c("G1[P1]"="goldenrod","G1[P2]"="orchid2"))
count_time_p2 <- sim$get_cells() %>% 
  dplyr::group_by(epistate) %>% 
  dplyr::summarise(n=n()) %>% 
  dplyr::mutate(prop = n / sum(n)) %>% 
  dplyr::mutate(time=sim$get_clock())

# Sample pre treatment
n_w <- n_h <- 25
ncells <- 0.8 * n_w * n_h
bbox <- sim$search_sample(c("G1" = ncells), n_w, n_h)
sim$sample_cells("S_PRE", bbox$lower_corner, bbox$upper_corner)
plot_tissue(sim,at_sample = "S_PRE",color_map = c("G1[P1]"="goldenrod","G1[P2]"="orchid2"))+facet_wrap(~epistate)



treatment_start <- sim$get_clock()
sim$set_rates(list("G1" = list(
  P1 = list(duplication = 0.5, death = 0.8, P2 = 0.2),
  P2 = list(duplication = 0.7, death = 0.3, P1 = 0.6)
)))

sim$run_up_to_time(100)

p3=plot_state(sim,color_map = c("G1[P1]"="goldenrod","G1[P2]"="orchid2"))
count_time_p3 <- sim$get_cells() %>% 
  dplyr::group_by(epistate) %>% 
  dplyr::summarise(n=n()) %>% 
  dplyr::mutate(prop = n / sum(n)) %>% 
  dplyr::mutate(time=sim$get_clock())



###############
sim$set_rates(list("G1" = list(
  P1 = list(duplication = 0.3, death = 0.9, P2 = 0.7),
  P2 = list(duplication = 0.7, death = 0.2, P1 = 0.1)
)))
sim$run_up_to_time(120)
# sim$get_counts()
p4=plot_state(sim,color_map = c("G1[P1]"="goldenrod","G1[P2]"="orchid2"))
count_time_p4 <- sim$get_cells() %>% 
  dplyr::group_by(epistate) %>% 
  dplyr::summarise(n=n()) %>% 
  dplyr::mutate(prop = n / sum(n)) %>% 
  dplyr::mutate(time=sim$get_clock())

###################
sim$set_rates(list("G1" = list(
  P1 = list(duplication = 0.9, death = 0.1, P2 = 0.4),
  P2 = list(duplication = 0.9, death = 0.1, P1 = 0.2)
)))
sim$run_up_to_time(140)

p5=plot_state(sim,color_map = c("G1[P1]"="goldenrod","G1[P2]"="orchid2"))
count_time_p5 <- sim$get_cells() %>% 
  dplyr::group_by(epistate) %>% 
  dplyr::summarise(n=n()) %>% 
  dplyr::mutate(prop = n / sum(n)) %>% 
  dplyr::mutate(time=sim$get_clock())

treatment_end <- sim$get_clock()
# Sample pre treatment
n_w <- n_h <- 25
ncells <- 0.8 * n_w * n_h
bbox <- sim$search_sample(c("G1" = ncells), n_w, n_h)
sim$sample_cells("S_POST_T1", bbox$lower_corner, bbox$upper_corner)

###################
sim$set_rates(list("G1" = list(
  P1 = list(duplication = 2, death = 0.1, P2 = 0.2), ### increase rates
  P2 = list(duplication = 2.3, death = 0.1, P1 = 0.5) ####### increase rates of growth like a relapse that growths faster
)))
sim$run_up_to_time(220)

p_ts <- plot_timeseries(sim,color_map = c("G1[P1]"="goldenrod","G1[P2]"="orchid2")) +
  annotate(
    "rect",
    xmin = treatment_start,
    xmax = treatment_end,
    ymin = -Inf,
    ymax = Inf,
    fill = "firebrick4",
    alpha = 0.2
  )
p7=plot_state(sim,color_map = c("G1[P1]"="goldenrod","G1[P2]"="orchid2"))
count_time_p7 <- sim$get_cells() %>% 
  dplyr::group_by(epistate) %>% 
  dplyr::summarise(n=n()) %>% 
  dplyr::mutate(prop = n / sum(n)) %>% 
  dplyr::mutate(time=sim$get_clock())
end_time_sim <- sim$get_clock()
start_time_sim <- 0
count_time_prop <- do.call("rbind",list(count_time_p1,count_time_p2,count_time_p3,count_time_p4,count_time_p6,count_time_p7))
pheno_prop_sim <- count_time_prop %>% 
  mutate(time=as.factor(round(time,0))) %>% 
  ggplot(
    aes(x = time, stratum = epistate, alluvium = epistate,
        y = prop, fill = epistate)) +
  geom_flow(alpha = 0.4) +
  geom_stratum(width = 1/3, color = "grey30") +
  geom_text(stat = "stratum", aes(label = epistate), size = 3) +
  annotate(
    "segment",
    x = 1,
    xend = 6,
    y = 1.02,
    yend = 1.02,
    linewidth = 10,
    colour = "firebrick4",
    alpha=0.2
  ) +
  scale_fill_manual(values=c("P1"="goldenrod","P2"="orchid2"))+
  theme_minimal()

# Sample pre treatment
n_w <- n_h <- 25
ncells <- 0.8 * n_w * n_h
bbox <- sim$search_sample(c("G1" = ncells), n_w, n_h)
sim$sample_cells("S_POST_T2", bbox$lower_corner, bbox$upper_corner)



###################
sim$set_rates(list("G1" = list(
  P1 = list(duplication = 1.5, death = 0.6, P2 = 0.4),
  P2 = list(duplication = 1.7, death = 0.4, P1 = 0.3))))

sim$run_up_to_time(240)
p8=plot_state(sim,color_map = c("G1[P1]"="goldenrod","G1[P2]"="orchid2"))
count_time_p8 <- sim$get_cells() %>% 
  dplyr::group_by(epistate) %>% 
  dplyr::summarise(n=n()) %>% 
  dplyr::mutate(prop = n / sum(n)) %>% 
  dplyr::mutate(time=sim$get_clock())


###################
sim$set_rates(list("G1" = list(
  P1 = list(duplication = 0.4, death = 0.6, P2 = 0.9),
  P2 = list(duplication = 0.6, death = 0.4, P1 = 0.1))))

sim$run_up_to_time(260)
p9=plot_state(sim,color_map = c("G1[P1]"="goldenrod","G1[P2]"="orchid2"))
count_time_p9 <- sim$get_cells() %>% 
  dplyr::group_by(epistate) %>% 
  dplyr::summarise(n=n()) %>% 
  dplyr::mutate(prop = n / sum(n)) %>% 
  dplyr::mutate(time=sim$get_clock())
plot_timeseries(simulation = sim)

