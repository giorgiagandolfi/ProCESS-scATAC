# Utility functions for the ProCESS single-treatment resistance simulation

get_state_counts <- function(sim, phase) {
  sim$get_cells() |>
    dplyr::count(epistate, name = "n") |>
    tidyr::complete(epistate = c("P1", "P2"), fill = list(n = 0L)) |>
    dplyr::mutate(
      total_cells = sum(n),
      prop = dplyr::if_else(total_cells > 0, n / total_cells, 0),
      time = as.numeric(sim$get_clock()),
      phase = phase
    )
}

run_and_record <- function(sim, end_time, phase, step = 2) {
  start_time <- as.numeric(sim$get_clock())
  if (end_time <= start_time) stop("end_time must be greater than current simulation time.")
  times <- seq(start_time + step, end_time, by = step)
  if (!length(times) || tail(times, 1) < end_time) times <- c(times, end_time)
  purrr::map_dfr(times, function(t) {
    sim$run_up_to_time(t)
    get_state_counts(sim, phase)
  })
}

sample_tissue_plot <- function(sim, sample_name, state_colors, n_w = 25, n_h = 25,
                               occupancy = 0.8) {
  ncells <- occupancy * n_w * n_h
  bbox <- sim$search_sample(c("G1" = ncells), n_w, n_h)
  sim$sample_cells(sample_name, bbox$lower_corner, bbox$upper_corner)
  plot_tissue(sim, at_sample = sample_name, color_map = state_colors) +
    ggplot2::facet_wrap(~epistate) +
    ggplot2::labs(title = sample_name)
}

add_treatment_window <- function(p, timeline, alpha = 0.08) {
  p +
    ggplot2::annotate("rect", xmin = timeline$treatment_start,
                      xmax = timeline$treatment_end,
                      ymin = -Inf, ymax = Inf, alpha = alpha) +
    ggplot2::geom_vline(xintercept = c(timeline$treatment_start,
                                       timeline$treatment_end),
                        linetype = "dashed", linewidth = 0.35)
}

prepare_treatment_dynamics <- function(count_time, timeline) {
  count_time |>
    dplyr::filter(phase == "Treatment") |>
    dplyr::mutate(treatment_time = time - timeline$treatment_start)
}

prepare_rate_history <- function(sim) {
  sim$get_rates_update_history() |>
    dplyr::mutate(
      time = as.numeric(time),
      epistate = factor(epistate, levels = c("P1", "P2")),
      transition = dplyr::case_when(
        event == "switch" & epistate == "P1" ~ "P1 → P2",
        event == "switch" & epistate == "P2" ~ "P2 → P1",
        TRUE ~ NA_character_
      )
    )
}

plot_total_population <- function(count_time, timeline) {
  dat <- count_time |>
    dplyr::select(time, phase, total_cells) |>
    dplyr::distinct() |>
    dplyr::arrange(time)
  p <- ggplot2::ggplot(dat, ggplot2::aes(time, total_cells)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 1.3) +
    ggplot2::labs(x = "Simulation time", y = "Total number of cells",
                  title = "Tumor burden during a single treatment") +
    ggplot2::theme_classic(base_size = 12)
  add_treatment_window(p, timeline)
}

plot_absolute_states <- function(count_time, timeline, epistate_colors) {
  p <- ggplot2::ggplot(count_time,
    ggplot2::aes(time, n, color = epistate, group = epistate)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 1.1) +
    ggplot2::scale_color_manual(values = epistate_colors,
      labels = c(P1 = "Sensitive (P1)", P2 = "Resistant (P2)")) +
    ggplot2::labs(x = "Simulation time", y = "Cell count", color = "Phenotype",
                  title = "Sensitive and resistant populations") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(legend.position = "top")
  add_treatment_window(p, timeline)
}

plot_composition <- function(count_time, timeline, epistate_colors) {
  p <- ggplot2::ggplot(count_time, ggplot2::aes(time, prop, fill = epistate)) +
    ggplot2::geom_area(position = "stack", alpha = 0.85, linewidth = 0.25) +
    ggplot2::scale_fill_manual(values = epistate_colors,
      labels = c(P1 = "Sensitive (P1)", P2 = "Resistant (P2)")) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1),
                                expand = c(0, 0)) +
    ggplot2::labs(x = "Simulation time", y = "Cell proportion", fill = "Phenotype",
                  title = "Phenotypic composition") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(legend.position = "top")
  add_treatment_window(p, timeline)
}

plot_treatment_burden <- function(treatment_dynamics, treatment_duration) {
  dat <- treatment_dynamics |>
    dplyr::select(treatment_time, total_cells) |>
    dplyr::distinct() |>
    dplyr::arrange(treatment_time) |>
    dplyr::mutate(relative_burden = total_cells / dplyr::first(total_cells))
  ggplot2::ggplot(dat, ggplot2::aes(treatment_time, relative_burden)) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dotted") +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::geom_point(size = 1.6) +
    ggplot2::scale_y_continuous(labels = scales::percent_format()) +
    ggplot2::scale_x_continuous(breaks = seq(0, treatment_duration, by = 5)) +
    ggplot2::labs(x = "Time since treatment start", y = "Relative tumor burden",
                  title = "Normalized treatment response",
                  subtitle = "Tumor burden at treatment start = 100%") +
    ggplot2::theme_classic(base_size = 12)
}

plot_p2_dynamics <- function(treatment_dynamics, treatment_duration) {
  dat <- treatment_dynamics |>
    dplyr::filter(epistate == "P2") |>
    dplyr::arrange(treatment_time) |>
    dplyr::mutate(delta_P2 = prop - dplyr::first(prop))
  p1 <- ggplot2::ggplot(dat, ggplot2::aes(treatment_time, prop)) +
    ggplot2::geom_line(linewidth = 1.1) + ggplot2::geom_point(size = 1.6) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
    ggplot2::scale_x_continuous(breaks = seq(0, treatment_duration, by = 5)) +
    ggplot2::labs(x = "Time since treatment start", y = "P2 proportion",
                  title = "Resistant phenotype during treatment") +
    ggplot2::theme_classic(base_size = 12)
  p2 <- ggplot2::ggplot(dat, ggplot2::aes(treatment_time, delta_P2)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted") +
    ggplot2::geom_line(linewidth = 1.1) + ggplot2::geom_point(size = 1.6) +
    ggplot2::scale_y_continuous(labels = scales::percent_format()) +
    ggplot2::scale_x_continuous(breaks = seq(0, treatment_duration, by = 5)) +
    ggplot2::labs(x = "Time since treatment start", y = "Change in P2 proportion",
                  title = "P2 enrichment relative to treatment start") +
    ggplot2::theme_classic(base_size = 12)
  list(proportion = p1, delta = p2)
}

plot_growth_rates <- function(rates_history, timeline, epistate_colors) {
  dat <- rates_history |> dplyr::filter(event %in% c("duplication", "death"))
  p <- ggplot2::ggplot(dat, ggplot2::aes(time, rate, color = epistate,
    linetype = event, group = interaction(epistate, event))) +
    ggplot2::geom_step(linewidth = 1.05, direction = "hv") +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::scale_color_manual(values = epistate_colors,
      labels = c(P1 = "P1 sensitive", P2 = "P2 resistant")) +
    ggplot2::labs(x = "Simulation time", y = "Event rate", color = "Phenotype",
                  linetype = "Event", title = "Duplication and death rates") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(legend.position = "top")
  add_treatment_window(p, timeline)
}

plot_switch_rates <- function(rates_history, timeline) {
  dat <- rates_history |> dplyr::filter(event == "switch")
  p <- ggplot2::ggplot(dat, ggplot2::aes(time, rate, color = transition,
                                         group = transition)) +
    ggplot2::geom_step(linewidth = 1.1, direction = "hv") +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::labs(x = "Simulation time", y = "Switching rate", color = "Transition",
                  title = "Phenotypic switching rates",
                  subtitle = "Single exposure: no treatment-history-dependent memory term") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(legend.position = "top")
  add_treatment_window(p, timeline)
}

save_plot <- function(plot, filename, width = 10, height = 6, dpi = 300) {
  ggplot2::ggsave(filename, plot = plot, width = width, height = height, dpi = dpi)
  invisible(filename)
}
