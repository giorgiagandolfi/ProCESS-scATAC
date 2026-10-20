library(ProCESS)
library(ggplot2)
library(dplyr)
library(patchwork)
library(ggalluvial)
library(tidyr)

setwd("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project/ProCESS-scATAC/scripts/simulate_fragments/peak_based/0_process_simulations/5_clones_evolution/")
source("../../utils.R")

sim <- recover_simulation("CRC_five_clones/")
sample_forest <- load_sample_forest("sample_forest.sff")



make_case_color_maps <- function() {
  list(
    epistate = c(
      "G1[P1]" = "goldenrod",
      "G2[P1]" = "goldenrod",
      "G3[P1]" = "goldenrod",
      "G4[P1]" = "goldenrod",
      "G5[P1]" = "goldenrod",
      
      "G1[P2]" = "orchid2",
      "G2[P2]" = "orchid2",
      "G3[P2]" = "orchid2",
      "G4[P2]" = "orchid2",
      "G5[P2]" = "orchid2",
      
      "G1[P3]" = "forestgreen",
      "G2[P3]" = "forestgreen",
      "G3[P3]" = "forestgreen",
      "G4[P3]" = "forestgreen",
      "G5[P3]" = "forestgreen"
    ),
    mutant = c(
      "G1[P1]"="coral2", "G2[P1]"="turquoise4", "G3[P1]"="darkorange","G4[P1]"="mediumorchid3","G5[P1]"="dodgerblue2",
      "G1[P2]"="coral2", "G2[P2]"="turquoise4", "G3[P2]"="darkorange","G4[P2]"="mediumorchid3","G5[P2]"="dodgerblue2",
      "G1[P3]"="coral2", "G2[P3]"="turquoise4", "G3[P3]"="darkorange","G4[P3]"="mediumorchid3","G5[P3]"="dodgerblue2"
    ),
    epi = c(P1="goldenrod", P2="orchid2", P3="forestgreen"),
    gen = c(G1="coral2", G2="turquoise4", G3="darkorange",G4="mediumorchid3",G5="dodgerblue2")
  )
}

color_map <- make_case_color_maps()
tissue_sample_1 <- plot_tissue(sim,before_sample = "S1",color_map = color_map$mutant)
tissue_sample_2 <- plot_tissue(sim,before_sample = "S2",color_map = color_map$mutant)
tissue_sample_3 <- plot_tissue(sim,before_sample = "S3",color_map = color_map$mutant)

plot_timeseries <- plot_timeseries(sim,color_map = color_map$mutant)
plot_forest_mutant <- plot_forest(sample_forest,color_map=color_map$mutant)
plot_forest_epistate <- plot_forest(sample_forest,color_map=color_map$epistate)

plot_epistate_alluvial <- function(df, colors) {
  ggplot2::ggplot(df,
                  ggplot2::aes(x = mutant, stratum = epistate, alluvium = epistate,
                               y = prop, fill = epistate)) +
    ggalluvial::geom_flow(alpha = 0.4) +
    ggalluvial::geom_stratum(width = 1/3, color = "grey30") +
    ggplot2::geom_text(stat = "stratum", ggplot2::aes(label = epistate), size = 3) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::scale_y_continuous(labels = scales::percent_format()) +
    ggplot2::labs(
                  x = "Genetic clone", y = "Within-clone proportion", fill = "Phenotype") +
    ggplot2::theme_minimal(base_size = 12)
}

plot_mutant_alluvial <- function(df, colors) {
  ggplot2::ggplot(df,
                  ggplot2::aes(x = epistate, stratum = mutant, alluvium = mutant,
                               y = prop, fill = mutant)) +
    ggalluvial::geom_flow(alpha = 0.4) +
    ggalluvial::geom_stratum(width = 1/3, color = "grey30") +
    ggplot2::geom_text(stat = "stratum", ggplot2::aes(label = mutant), size = 3) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::scale_y_continuous(labels = scales::percent_format()) +
    ggplot2::labs(
                  x = "Phenotype", y = "Within-phenotype proportion", fill = "Genetic clone") +
    ggplot2::theme_minimal(base_size = 12)
}
summarise_sample <- function(sample_forest) {
  nodes <- sample_forest$get_nodes() %>% dplyr::filter(!is.na(sample))
  
  epi_by_mutant <- nodes %>%
    dplyr::count(mutant, epistate, name = "n") %>%
    dplyr::group_by(mutant) %>%
    dplyr::mutate(prop = n / sum(n)) %>%
    dplyr::ungroup()
  
  mutant_by_epi <- nodes %>%
    dplyr::count(epistate, mutant, name = "n") %>%
    dplyr::group_by(epistate) %>%
    dplyr::mutate(prop = n / sum(n)) %>%
    dplyr::ungroup()
  
  list(epi_by_mutant = epi_by_mutant, mutant_by_epi = mutant_by_epi)
}
sample_summary <- summarise_sample(sample_forest)
sankey_epi <- plot_epistate_alluvial(sample_summary$epi_by_mutant, color_map$epi)
sankey_gen <- plot_mutant_alluvial(sample_summary$mutant_by_epi, color_map$gen)
sample_prop <- sample_forest$get_nodes() %>% 
  dplyr::filter(!is.na(sample)) %>% 
  group_by(sample) %>% 
  count(mutant,epistate)
plot_df <- sample_prop %>%
  pivot_longer(
    cols = c(mutant, epistate),
    names_to = "composition",
    values_to = "category"
  ) %>%
  group_by(sample, composition, category) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  group_by(sample, composition) %>%
  mutate(
    fraction = n / sum(n),
    label = if_else(
      fraction >= 0.03,
      paste0(
        category, "\n",
        scales::percent(fraction, accuracy = 0.1)
      ),
      ""
    )
  ) %>%
  ungroup() %>%
  mutate(
    composition = factor(
      composition,
      levels = c("mutant", "epistate"),
      labels = c("Genetic clone", "Phenotype")
    ),
    category = factor(
      category,
      levels = c(paste0("G", 1:5), paste0("P", 1:3))
    )
  )

sample_composition_plots <- ggplot(plot_df, aes(x = "", y = fraction, fill = category)) +
  geom_col(width = 1, colour = "white", linewidth = 0.4) +
  geom_text(
    aes(label = category),
    position = position_stack(vjust = 0.5),
    size = 3.5
  ) +
  coord_polar(theta = "y") +
  facet_grid(sample~composition) +
  scale_fill_manual(
    values = c(
      color_map$gen,
      color_map$epi
    )
  ) +
  labs(fill = "Clone / phenotype") +
  theme_void()+
  theme(
    strip.text = element_text(size = 12),
    legend.position = "bottom"
  )

transition_graph_g1 <- plot_mutant_transition_graph(sim = sim,mutant_name = "G1",colors = color_map$epi)
transition_graph_g2 <-plot_mutant_transition_graph(sim = sim,mutant_name = "G2",colors = color_map$epi)
transition_graph_g3 <-plot_mutant_transition_graph(sim = sim,mutant_name = "G3",colors = color_map$epi)
transition_graph_g4 <-plot_mutant_transition_graph(sim = sim,mutant_name = "G4",colors = color_map$epi)
transition_graph_g5  <-plot_mutant_transition_graph(sim = sim,mutant_name = "G5",colors = color_map$epi)


plot_tissue_evolution <- wrap_plots(list(plot_timeseries,tissue_sample_1,tissue_sample_2,tissue_sample_3,
                                         sample_composition_plots),
                # transition_graph_g1,transition_graph_g2,transition_graph_g3,transition_graph_g4,transition_graph_g5),
           design = "aaaee\nbcdee",guides="collect")+plot_annotation(tag_levels = "A") & theme(legend.position = "none")

ggsave(filename = "plot_tissue_evolution_case5.png",plot = plot_tissue_evolution,width = 9,height = 6,dpi = 300)

wrap_plots(list(plot_forest_mutant,plot_forest_epistate,sankey_epi,sankey_gen),
           design = "aabb\naabb\naabb\nddee",guides="collect") +plot_annotation(tag_levels = "A") & theme(legend.position = "none")
ggsave(filename = "plot_forest_case5.png",width = 9,height = 6,dpi = 300)
