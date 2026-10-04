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
                        width = 2000, height = 2000,save_directory = T,name = "case3")
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
  P1 = list(duplication = 0.4, death = 0.2, P2 = 0.8),
  P2 = list(duplication = 0.8, death = 0.1, P1 = 0.1)
)))
sim$run_up_to_time(140)
p5=plot_state(sim,color_map = c("G1[P1]"="goldenrod","G1[P2]"="orchid2"))
count_time_p5 <- sim$get_cells() %>% 
  dplyr::group_by(epistate) %>% 
  dplyr::summarise(n=n()) %>% 
  dplyr::mutate(prop = n / sum(n)) %>% 
  dplyr::mutate(time=sim$get_clock())

###################
sim$set_rates(list("G1" = list(
  P1 = list(duplication = 0.9, death = 0.1, P2 = 0.4),
  P2 = list(duplication = 0.9, death = 0.1, P1 = 0.2)
)))
sim$run_up_to_time(160)

p6=plot_state(sim,color_map = c("G1[P1]"="goldenrod","G1[P2]"="orchid2"))
count_time_p6 <- sim$get_cells() %>% 
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
  P1 = list(duplication = 1.7, death = 0.1, P2 = 0.2), ### increase rates
  P2 = list(duplication = 1.5, death = 0.1, P1 = 0.5) ####### increase rates of growth like a relapse that growths faster
)))
sim$run_up_to_time(200)

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
  # labs(title = "Epistate proportions across mutant status",
  #      x = "Mutant", y = "Count", fill = "Epistate") +
  scale_fill_manual(values=c("P1"="goldenrod","P2"="orchid2"))+
  theme_minimal()

# Sample pre treatment
n_w <- n_h <- 25
ncells <- 0.8 * n_w * n_h
bbox <- sim$search_sample(c("G1" = ncells), n_w, n_h)
sim$sample_cells("S_POST_T2", bbox$lower_corner, bbox$upper_corner)


# df_treatment <- data.frame(
#   start = c(start_time_sim, treatment_start, treatment_end),
#   end   = c(treatment_start, treatment_end, end_time_sim),
#   period = c("No treatment", "Treatment", "No treatment")
# )
# tratment_timeline_plt <-ggplot(df_treatment) +
#   geom_rect(
#     aes(
#       xmin = start,
#       xmax = end,
#       ymin = 0,
#       ymax = 1,
#       fill = period
#     )
#   ) +
#   scale_fill_manual(
#     values = c(
#       "No treatment" = "gainsboro",
#       "Treatment" = "firebrick3"
#     )
#   ) +
#   scale_x_continuous(
#     limits = c(start_time_sim, end_time_sim),
#     breaks = seq(start_time_sim, end_time_sim, by = 20)
#   ) +
#   labs(
#     x = "Simulation time",
#     y = NULL,
#     fill = NULL
#   ) +
#   theme_minimal() +
#   theme(
#     axis.text.y = element_blank(),
#     axis.ticks.y = element_blank()
#   )


#sim <- recover_simulation("case3/")
sample_forest <- sim$get_sample_forest()
sample_forest$save("sample_forest_atac_case_3.sff")
plot_sample_forest <- plot_forest(sample_forest,color_map = c("G1[P1]"="goldenrod","G1[P2]"="orchid2"))+
  annotate(
    "rect",
    ymin = end_time_sim-treatment_end,
    ymax = treatment_start,
    xmin = -Inf,
    xmax = Inf,
    fill = "firebrick4",
    alpha = 0.2
  )
f_plot1 <-wrap_plots(list(p1,p2,p3,p4,p5,p6,p7,
                         p_ts,pheno_prop_sim),design = "abcdefg\nhhhhhhh\niiiiiii",guides = "collect")&theme(legend.position = "bottom")



rates_history <- sim$get_rates_update_history()
library(ggplot2)
library(dplyr)

switch_rates <- rates_history %>%
  filter(event == "switch") %>%
  mutate(
    transition = case_when(
      epistate == "P1" ~ "P1 \u2192 P2",
      epistate == "P2" ~ "P2 \u2192 P1"
    ),
    time_label = factor(
      round(time),
      levels = c(0, 80, 100, 120, 140, 160)
    )
  )

switch_rates_plt <- ggplot(
  switch_rates,
  aes(x = time_label, y = transition, fill = rate)
) +
  
  # Heatmap
  geom_tile(
    color = "white",
    linewidth = 1.5
  ) +
  
  # Switching-rate values
  geom_text(
    aes(label = sprintf("%.2f", rate)),
    size = 5,
    # fontface = "bold"
  ) +
  
  # Treatment bar
  annotate(
    "segment",
    x = 2,
    xend = 6,
    y = 2.65,
    yend = 2.65,
    linewidth = 5,
    colour = "firebrick4",
    alpha=0.2
  ) +
  
  scale_fill_gradient(
    low = "white",
    high = "steelblue4",
    limits = c(0, 1),
    name = "Switching rate"
  ) +
  
  scale_y_discrete(
    expand = expansion(mult = c(0.1, 0.45))
  ) +
  
  labs(
    x = "Simulation time",
    y = NULL
  ) +
  
  coord_cartesian(clip = "off") +
  
  theme_minimal() +
  
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(
      size = 13,
      # face = "bold"
    ),
    legend.position = "right",
    plot.margin = margin(30, 10, 10, 10)
  )

f_plot2 <-wrap_plots(list(p_ts,pheno_prop_sim,switch_rates_plt,plot_sample_forest),
                     design = "aaddd\nbbddd\nccddd",guides = "collect")&theme(legend.position = "bottom")



##################### MUTATION PART ##############
setwd("process_references_v1.3.5")

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
m_engine$change_rates_from(treatment_end+10, "G1", list("P2" = c(SNV = 5e-9, indel = 2e-9)))
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
  time = treatment_start,
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
# POST-TREATMENT
# chemotherapy-associated mutational processes switched off
# ============================================================

m_engine$add_exposure(
  time = treatment_end,
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
phylo_forest <- m_engine$place_mutations(sample_forest, num_of_preneoplatic_SNVs=800, num_of_preneoplatic_indels=200)

phylo_forest$save("../phylo_forest_atac_case_3.sff")




exp_df <- phylo_forest$get_exposures() %>%
  mutate(
    type = recode(
      type,
      "SNV"   = "SBS",
      "indel" = "ID"
    )
  )
# Repeat the final exposure state at the end of simulation
last_exposure <- exp_df %>%
  group_by(signature, type) %>%
  slice_max(time, n = 1, with_ties = FALSE) %>%
  mutate(time = 200)

exp_plot <- bind_rows(
  exp_df,
  last_exposure
)
exposure_time_plt <-ggplot(
  exp_plot,
  aes(
    x = time,
    y = exposure,
    color = signature,
    group = signature
  )
) +
  
  annotate(
    "rect",
    xmin = treatment_start,
    xmax = treatment_end,
    ymin = -Inf,
    ymax = Inf,
    fill = "firebrick4",
    alpha = 0.07
  ) +
  
  geom_step(
    direction = "hv",
    linewidth = 1.2
  ) +
  
  geom_point(
    data = exp_df,
    size = 2.8
  ) +
  
  facet_wrap(
    ~ type,
    ncol = 1
  ) +
  
  scale_x_continuous(
    breaks = seq(0, 200, 20),
    limits = c(0, 200),
    expand = c(0, 0)
  ) +
  
  scale_y_continuous(
    limits = c(0, 0.45),
    breaks = seq(0, 0.4, 0.1)
  ) +
  
  labs(
    x = "Simulation time",
    y = "Signature exposure",
    color = "Signature"
  ) +
  
  theme_minimal() +
  
  theme(
    legend.position = "right"
  )

species_info <- m_engine$get_species_info()

# Convert SNV and indel rates to long format
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
mut_rate_hist_plt <-ggplot(
  mut_rate_df,
  aes(
    x = time,
    y = mutation_rate_log,
    color = epistate,
    group = epistate
  )
) +
  
  geom_step(
    direction = "hv",
    linewidth = 1.2
  ) +
  
  geom_point(
    size = 3
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
  
  scale_x_continuous(
    breaks = seq(0, 200, 20),
    limits = c(0, 200)
  ) +
  annotate(
    "rect",
    xmin = treatment_start,
    xmax = treatment_end,
    ymin = -Inf,
    ymax = Inf,
    fill = "firebrick4",
    alpha = 0.07
  )+
  labs(
    x = "Simulation time",
    y = "Mutation rate",
    color = "Epistate"
  ) +
  theme_minimal()+
  
  theme(
    legend.position = "top",
  )


f_plot2 <-wrap_plots(list(p_ts,pheno_prop_sim,switch_rates_plt,exposure_time_plt,mut_rate_hist_plt,plot_sample_forest),
                     design = "aafff\nbbfff\nccfff\nddfff\neefff",guides = "collect")&theme(legend.position = "bottom")
f_plot2
