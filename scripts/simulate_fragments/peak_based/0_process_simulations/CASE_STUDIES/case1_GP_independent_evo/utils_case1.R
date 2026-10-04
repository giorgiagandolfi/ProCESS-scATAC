# Utilities for ProCESS case 2: genetic + epigenetic evolution

make_case1_color_maps <- function() {
  list(
    epistate = c(
      "G1[P1]" = "goldenrod",
      "G2[P1]" = "goldenrod",
      "G3[P1]" = "goldenrod",
      
      "G1[P2]" = "orchid2",
      "G2[P2]" = "orchid2",
      "G3[P2]" = "orchid2",
      
      "G1[P3]" = "forestgreen",
      "G2[P3]" = "forestgreen",
      "G3[P3]" = "forestgreen"
    ),
    mutant = c(
      "G1[P1]"="coral2", "G2[P1]"="turquoise4", "G3[P1]"="darkorange",
      "G1[P2]"="coral2", "G2[P2]"="turquoise4", "G3[P2]"="darkorange",
      "G1[P3]"="coral2", "G2[P3]"="turquoise4", "G3[P3]"="darkorange"
    ),
    epi = c(P1="goldenrod", P2="orchid2", P3="forestgreen"),
    gen = c(G1="coral2", G2="turquoise4", G3="darkorange")
  )
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

plot_epistate_alluvial <- function(df, colors) {
  ggplot2::ggplot(df,
    ggplot2::aes(x = mutant, stratum = epistate, alluvium = epistate,
                 y = prop, fill = epistate)) +
    ggalluvial::geom_flow(alpha = 0.4) +
    ggalluvial::geom_stratum(width = 1/3, color = "grey30") +
    ggplot2::geom_text(stat = "stratum", ggplot2::aes(label = epistate), size = 3) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::scale_y_continuous(labels = scales::percent_format()) +
    ggplot2::labs(title = "Phenotype composition by genetic clone",
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
    ggplot2::labs(title = "Genetic-clone composition by phenotype",
                  x = "Phenotype", y = "Within-phenotype proportion", fill = "Clone") +
    ggplot2::theme_minimal(base_size = 12)
}

plot_sample_box <- function(sim, lower, width = 40) {
  upper <- lower + width
  plot_tissue(sim) +
    ggplot2::geom_rect(xmin = lower[1], xmax = upper[1],
                       ymin = lower[2], ymax = upper[2],
                       fill = NA, color = "black", linewidth = 0.8) +
    ggplot2::facet_wrap(~mutant)
}

capture_q_plot <- function(q_obj) {
  grid::grid.grabExpr(ComplexHeatmap::draw(q_obj$plot_q))
}

make_case1_report_plot <- function(q_G1, q_G2, q_G3, sankey_gen, sankey_epi,
                                   p_forest_gen, p_forest_epi) {
  h1 <- capture_q_plot(q_G1); h2 <- capture_q_plot(q_G2); h3 <- capture_q_plot(q_G3)
  patchwork::wrap_plots(
    list(h1, h2, h3, sankey_gen, sankey_epi, p_forest_gen, p_forest_epi),
    design = "AABBCC\nDDDEEE\nFFFGGG\nFFFGGG"
  ) + patchwork::plot_annotation(tag_levels = "a") &
    ggplot2::theme(legend.position = "bottom")
}
