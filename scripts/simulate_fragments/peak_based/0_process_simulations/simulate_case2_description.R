rm(list=ls())
library(ProCESS)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggalluvial)
library(patchwork)
library(ComplexHeatmap)
setwd("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project/ProCESS-scATAC/scripts/simulate_fragments/peak_based/0_process_simulations")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project/ProCESS-scATAC/scripts/simulate_fragments/peak_based/utils.R")
# set.seed(0)
set.seed(999)
# sim <- TissueSimulation(epigenetic_states = c("E1", "E2","E3"),
#                         width = 1000, height = 1000,save_snapshots = F,name = "case2")


sim <- TissueSimulation(epigenetic_states = c("E1", "E2"),
                        width = 1000, height = 1000,save_snapshots = F,name = "case2")

# rates_genetic_clone_0 <- list(
#   E1 = list(duplication = 0.8, death = 0.1, E2 = 0.2,E3 = 0.6),
#   E2 = list(duplication = 0.8, death = 0.1, E1 = 0.3,E3=0.6),
#   E3 = list(duplication = 0.8, death = 0.1, E1 = 0.3,E2 = 0.2)
# )
# 
# 
# # add a mutant "A" and set its species rates
# rates_genetic_clone_A <- list(
#   E1 = list(duplication = 1, death = 0.1,E2 = 0.2,E3 = 0.6),
#   E2 = list(duplication = 1, death = 0.1, E1 = 0.3,E3=0.6),
#   E3 = list(duplication = 1, death = 0.1, E1 = 0.3,E2=0.2)
# )

# rates_genetic_clone_0 <- list(
#   E1 = list(duplication = 0.8, death = 0.1, E2 = 0.2,E3 = 0.6),
#   E2 = list(duplication = 0.8, death = 0.1, E1 = 0.3,E3=0.6),
#   E3 = list(duplication = 0.8, death = 0.1, E1 = 0.3,E2 = 0.2)
# )
# 
# 
# # add a mutant "A" and set its species rates
# rates_genetic_clone_A <- list(
#   E1 = list(duplication = 1, death = 0.1,E2 = 0.9,E3 = 0.1),
#   E2 = list(duplication = 1, death = 0.1, E1 = 0.05,E3=0.2),
#   E3 = list(duplication = 1, death = 0.1, E1 = 0.4,E2=0.8)
# )



### two phenotypes only
rates_genetic_clone_0 <- list(
  E1 = list(duplication = 0.8, death = 0.1, E2 = 0.1),
  E2 = list(duplication = 0.8, death = 0.1, E1 = 0.5)
  # E3 = list(duplication = 0.8, death = 0.1, E1 = 0.3,E2 = 0.2)
)


# add a mutant "A" and set its species rates
rates_genetic_clone_A <- list(
  E1 = list(duplication = 1, death = 0.1,E2 = 0.9),
  E2 = list(duplication = 1, death = 0.1, E1 = 0.05)
  # E3 = list(duplication = 1, death = 0.1, E1 = 0.4,E2=0.8)
)






# q_G1=get_epigenetic_Q(epigenetic_rates = rates_genetic_clone_0,epistate_name = "G1")
# q_G2=get_epigenetic_Q(epigenetic_rates = rates_genetic_clone_A,epistate_name = "G2")
# q_G1$plot_q
# q_G2$plot_q
# q_G1$steady_states
# q_G2$steady_states
# transition_matrix_plot <-q_G1$plot_q
sim$add_mutant("G1", rates_genetic_clone_0)
sim$place_cell("G1[E1]", 500, 500)

# let the simulation evolve until the species "A[E2]" has less than 10 cells
sim$run_up_to_time(20)


sim$add_mutant("G2",rates_genetic_clone_A)
starting_cell=sim$choose_border_cell_in("G1")
plot_tissue(sim)+
  geom_point(data = starting_cell,aes(x=position_x,y=position_y),inherit.aes = FALSE)+
  facet_wrap(~epistate)

sim$mutate_progeny(starting_cell,"G2")
sim$set_rates(list("G1" = list(duplication = 0.6)))
sim$run_up_to_time(60)
plot_tissue(sim)+
  facet_wrap(~epistate)




bbox_width=20
bbox1_p <- c(400,410)
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
  )+
  facet_wrap(~mutant)


sim$sample_cells("S3", bbox1_p, bbox1_q)
sample_forest <- sim$get_sample_forest()


sampled_epigenetic_prop= sample_forest$get_nodes() %>%
  filter(!is.na(sample)) %>%
  group_by(epistate,mutant) %>%
  summarise(n=n())
sampled_epigenetic_prop <- sampled_epigenetic_prop %>%
  group_by(mutant) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

sankey_epi <-ggplot(sampled_epigenetic_prop,
                    aes(x = mutant, stratum = epistate, alluvium = epistate,
                        y = prop, fill = epistate)) +
  geom_flow(alpha = 0.4) +
  geom_stratum(width = 1/3, color = "grey30") +
  geom_text(stat = "stratum", aes(label = epistate), size = 3) +
  labs(title = "Epistate proportions across mutant status",
       x = "Mutant", y = "Count", fill = "Epistate") +
  scale_fill_manual(values=c("E1"="forestgreen","E2"="goldenrod","E3"="orchid2"))+
  theme_minimal()
sampled_genetic_prop= sample_forest$get_nodes() %>%
  filter(!is.na(sample)) %>%
  group_by(epistate,mutant) %>%
  summarise(n=n())
sampled_genetic_prop <- sampled_genetic_prop %>%
  group_by(epistate) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

sankey_gen <-ggplot(sampled_genetic_prop,
                    aes(x = epistate, stratum = mutant, alluvium = mutant,
                        y = prop, fill = mutant)) +
  geom_flow(alpha = 0.4) +
  geom_stratum(width = 1/3, color = "grey30") +
  geom_text(stat = "stratum", aes(label = mutant), size = 3) +
  labs(title = "Mutant proportions across epistate",
       x = "Epistate", y = "Count", fill = "Mutant") +
  scale_fill_manual(values=c("G1"="coral2","G2"="turquoise4","G3"="darkorange"))+
  theme_minimal()
sankey_gen+sankey_epi

color_map_epi = c("G1[E1]"="forestgreen",
                  "G2[E1]"="forestgreen",
                  "G3[E1]"="forestgreen",
                  "G1[E2]"="goldenrod",
                  "G2[E2]"="goldenrod",
                  "G3[E2]"="goldenrod",
                  "G1[E3]"="orchid2",
                  "G2[E3]"="orchid2",
                  "G3[E3]"="orchid2"
)
color_map_gen = c("G1[E1]"="coral2",
                  "G2[E1]"="turquoise4",
                  "G3[E1]"="darkorange",
                  "G1[E2]"="coral2",
                  "G2[E2]"="turquoise4",
                  "G3[E2]"="darkorange",
                  "G1[E3]"="coral2",
                  "G2[E3]"="turquoise4",
                  "G3[E3]"="darkorange"
)
p1_forest<-plot_forest(sample_forest,color_map = color_map_epi)+
  theme(legend.position = "none")+ggtitle(label="Cell coloured by epistate")
p2_forest <-plot_forest(sample_forest,color_map = color_map_gen)+
  theme(legend.position = "none")+ggtitle(label="Cell coloured by genetic clone")
p2_forest+p1_forest





plot_data <- sampled_epigenetic_prop %>%
  group_by(mutant) %>%
  mutate(
    percent = prop / sum(prop) * 100,
    label = paste0(
      epistate,
      "\n",
      round(percent, 1),
      "%"
    )
  ) %>%
  ungroup()

ggplot(
  plot_data,
  aes(
    x = "",
    y = prop,
    fill = epistate
  )
) +
  geom_col(
    width = 1,
    color = "white"
  ) +
  coord_polar(theta = "y") +
  facet_wrap(~ mutant) +
  scale_fill_manual(values=c("E1"="forestgreen","E2"="goldenrod","E3"="orchid2"))+
  # geom_text(
  #   aes(
  #     label = label
  #   ),
  #   position = position_stack(
  #     vjust = 0.5
  #   ),
  #   size = 4
  # ) +
  # labs(
  #   title = "Epigenetic state composition by mutant group",
  #   fill = "Epigenetic state"
  # ) +
  theme_void()
