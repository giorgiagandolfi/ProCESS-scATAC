library(dplyr)
library(ggplot2)

prj_dir <- "/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project/ProCESS-scATAC/scripts/simulate_fragments/peak_based/0_process_simulations/CASE_STUDIES/case4_epigenetic_memory/"
output_dir <- file.path(
  prj_dir,
  "comparison_memory"
)

memory <- read.csv(
  file.path(
    output_dir,
    "memory_summary_stats.csv"
  )
)

no_memory <- read.csv(
  file.path(
    output_dir,
    "no_memory_summary_stats.csv"
  )
)


comparison <- bind_rows(
  no_memory,
  memory
)


# =============================================================================
# PERCENT REDUCTION IN RELAPSE TIME
# =============================================================================

t_no_memory <- comparison %>%
  filter(model == "No memory") %>%
  pull(time_to_relapse_T2)

t_memory <- comparison %>%
  filter(model == "Epigenetic memory") %>%
  pull(time_to_relapse_T2)


percent_reduction <- 100 * (
  t_no_memory - t_memory
) / t_no_memory


cat(
  "Reduction in T2 relapse time:",
  round(percent_reduction, 1),
  "%\n"
)

comparison %>%
  ggplot(
    aes(
      x = model,
      y = time_to_relapse_T2,
      fill = model
    )
  ) +
  
  geom_col(
    width = 0.4
  ) +
  
  geom_text(
    aes(
      label = sprintf(
        "%.2f",
        time_to_relapse_T2
      )
    ),
    vjust = -0.5,
    size = 5
  ) +
  
  annotate(
    "text",
    x = 1.5,
    y = max(comparison$time_to_relapse_T2) * 1.15,
    label = sprintf(
      "%.1f%% reduction",
      percent_reduction
    ),
    size = 5,
    # fontface = "bold"
  ) +
  
  scale_fill_manual(
    values = c(
      "No memory" = "grey70",
      "Epigenetic memory" = "firebrick2"
    )
  ) +
  
  
  labs(
    x = NULL,
    y = "Time to resistant relapse",
    title = "Effect of epigenetic memory on rechallenge"
  ) +
  
  theme_minimal(
    base_size = 14
  ) +
  
  theme(
    legend.position = "none",
    plot.title = element_text(
      hjust = 0.5,
      # face = "bold"
    )
  )

plot_timeseries_smooth <- function(
    simulation,
    color_map = NULL,
    span = 0.3,
    show_points = TRUE
) {
  
  stopifnot(
    inherits(simulation, "Rcpp_TissueSimulation")
  )
  
  counts <- simulation$get_count_history() %>%
    ProCESS:::add_species_col()
  
  time <- simulation$get_clock() %>%
    round(digits = 3)
  
  if (is.null(color_map)) {
    color_map <- ProCESS:::get_species_colors(simulation)
  }
  
  counts$species <- factor(
    counts$species,
    levels = unique(names(color_map))
  )
  
  p <- ggplot(
    counts,
    aes(
      x = time,
      y = count,
      color = species
    )
  ) +
    
    # Original trajectory, lightly shown
    geom_line(
      alpha = 0.25,
      linewidth = 0.5
    ) +
    
    # Smooth trajectory
    geom_smooth(
      method = "loess",
      formula = y ~ x,
      se = FALSE,
      span = span,
      linewidth = 1.2
    ) +
    
    labs(
      x = "Time",
      y = "Cell count",
      color = "Species"
    ) +
    
    ProCESS:::my_theme() +
    
    scale_color_manual(
      values = color_map
    ) +
    
    theme(
      legend.position = "bottom"
    )
  
  if (show_points) {
    
    p <- p +
      geom_point(
        alpha = 0.25,
        size = 0.8
      )
  }
  
  return(p)
}


plot_timeseries_start0 <- function(
    simulation,
    color_map = NULL
) {
  counts <- simulation$get_count_history() %>% ProCESS:::add_species_col()
  counts <- counts %>%
    bind_rows(
      counts %>%
        distinct(mutant, epistate, species) %>%
        mutate(
          count = 0L,
          time = 0
        )
    ) %>%
    arrange(time, mutant, epistate)
  time <- simulation$get_clock() %>% round(digits = 3)
  if (is.null(color_map)) {
    color_map <- get_species_colors(simulation)
  }
  counts$species <- factor(counts$species, levels = unique(names(color_map)))
  ggplot2::ggplot(counts) + ggplot2::geom_line(ggplot2::aes(x = time, 
                                                            y = .data$count, color = .data$species)) + ggplot2::geom_point(ggplot2::aes(x = time, 
                                                                                                                                        y = .data$count, color = .data$species)) + ggplot2::labs(x = "Time", 
                                                                                                                                                                                                 y = "Cell count", color = "Species") + my_theme() + ggplot2::scale_color_manual(values = color_map) + 
    ggplot2::theme(legend.position = "bottom")
}
plot_firings_my <- function(simulation){
  {
    firings <- simulation$get_firings() %>% ProCESS:::add_species_col()
    plot <- ggplot2::ggplot(firings) + ggplot2::geom_bar(stat = "identity", 
                                                         ggplot2::aes(x = "", y = .data$fired, fill = .data$event))
    if ("epistate" %in% colnames(firings)) {
      plot <- plot + ggplot2::facet_grid(.data$mutant ~ .data$epistate)
    }
    else {
      plot <- plot + ggplot2::facet_grid(.data$mutant)
    }
    plot + ggplot2::labs(fill = "Event",x = NULL, y = NULL) + ProCESS:::my_theme() + ggplot2::scale_fill_brewer(palette = "Dark2") + 
      ggplot2::theme(legend.position = "bottom")
  }
}

