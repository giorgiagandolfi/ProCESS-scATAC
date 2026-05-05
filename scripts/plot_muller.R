
plot_Muller <- function(simulation, color_map = NULL) {
  
  stopifnot(inherits(simulation, "Rcpp_TissueSimulation"))
  
  # -------------------------
  # Population history
  # -------------------------
  df_populations <- simulation$get_count_history() %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(Identity = paste0(.data$mutant, .data$epistate)) %>%
    dplyr::rename(
      Generation = .data$time,
      Population = .data$count
    ) %>%
    dplyr::select(.data$Generation, .data$Identity, .data$Population)
  
  # -------------------------
  # Lineage edges
  # -------------------------
  df_edges <- simulation$get_lineage_graph() %>%
    dplyr::distinct(.data$ancestor, .data$progeny) %>%
    dplyr::rename(
      Parent   = .data$ancestor,
      Identity = .data$progeny
    ) %>%
    dplyr::select(.data$Parent, .data$Identity)
  
  df_edges <- ProCESS:::collapse_loops(df_edges)
  
  # -------------------------
  # Wild-type padding
  # -------------------------
  max_tumour_size <- df_populations %>%
    dplyr::group_by(.data$Generation) %>%
    dplyr::summarise(Population = sum(.data$Population), .groups = "drop") %>%
    dplyr::pull(.data$Population) %>%
    max()
  
  max_tumour_size <- max_tumour_size * 1.05
  
  wt_dynamics <- df_populations %>%
    dplyr::group_by(.data$Generation) %>%
    dplyr::summarise(Population = sum(.data$Population), .groups = "drop") %>%
    dplyr::mutate(
      Identity   = "Wild-type",
      Population = max_tumour_size - .data$Population
    )
  
  t_wt_dynamics <- dplyr::bind_rows(wt_dynamics, df_populations)
  
  # -------------------------
  # ✅ Clone ordering by first appearance
  # -------------------------
  clone_order <- df_populations %>%
    dplyr::filter(.data$Population > 0) %>%
    dplyr::group_by(.data$Identity) %>%
    dplyr::summarise(
      first_time = min(.data$Generation),
      .groups = "drop"
    ) %>%
    dplyr::arrange(.data$first_time) %>%
    dplyr::pull(.data$Identity)
  
  t_wt_dynamics <- t_wt_dynamics %>%
    dplyr::mutate(
      Identity = factor(
        .data$Identity,
        levels = c("Wild-type", clone_order)
      )
    )
  
  # -------------------------
  # Colors
  # -------------------------
  if (is.null(color_map)) {
    color_map <- ProCESS:::get_species_colors(simulation$get_species())
  }
  
  group_name <- ProCESS:::get_group_cell_name(simulation)
  
  # -------------------------
  # Muller plot
  # -------------------------
  suppressWarnings({
    muller_df <- ggmuller::get_Muller_df(df_edges, t_wt_dynamics)
    
    plot <- ggmuller::Muller_pop_plot(
      muller_df,
      add_legend = TRUE,
      palette = c(`Wild-type` = "gainsboro", color_map)
    ) +
      ProCESS:::my_theme() +
      ggplot2::guides(fill = ggplot2::guide_legend(group_name))
  })
  
  return(plot)
}


## This file is part of the ProCESS (https://github.com/caravagnalab/ProCESS/).
## Copyright (C) 2023-2025 - Giulio Caravagna <gcaravagna@units.it>
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <https://www.gnu.org/licenses/>.

#' Annotate a plot of cell divisions
#'
#' @description
#' It annotates a plot of cell divisions with information from sampling
#' times and MRCAs for all available samples
#'
#' @param tree_plot The output of `plot_forest`.
#' @param forest The original forest object from which the input to
#'        `plot_forest`
#' has been derived.
#' @param samples If `TRUE` it annotates samples.
#' @param MRCAs If `TRUE` it annotates MRCAs.
#' @param exposures If `TRUE` it annotates exposures to mutational signatures.
#' @param facet_signatures If `TRUE` and if `exposures` is `TRUE` it creates a
#'   faceted forest plot where the exposure to each signature is annotated on
#'   a separated plot.
#' @param drivers If `TRUE` it annotates drivers on the node they originated.
#' @param add_driver_label If `TRUE` and if `drivers` is `TRUE` it annotates
#'   the driver name.
#'
#' @return A `ggraph` tree plot.
#' @export
#'
#' @examples
#' sim <- TissueSimulation()
#' sim$add_mutant(name = "A", growth_rates = 0.08, death_rates = 0.01)
#' sim$place_cell("A", 500, 500)
#' sim$run_up_to_time(60)
#' sim$sample_cells("MySample", c(500, 500), c(510, 510))
#' m_engine = MutationEngine(setup_code = "demo")
#'
#' m_engine$add_mutant(mutant_name = "A",
#'                     passenger_rates = c(SNV = 1e-9),
#'                     drivers = list(SNV("22", 10510210, "C"),
#'                                    CNA(type = "A", "22", chr_pos = 10303470,
#'                                        len = 200000)))
#' m_engine$add_exposure(coefficients = c(SBS13 = 0.2, SBS1 = 0.8))
#' m_engine$add_exposure(time=50, coefficients = c(SBS17b = 0.2, SBS3 = 0.8))
#'
#' forest = sim$get_sample_forest()
#' forest$get_samples_info()
#' forest_muts = m_engine$place_mutations(forest, 1000, 500)
#' tree_plot = plot_forest(forest)
#' annotate_forest(tree_plot, forest_muts, samples = T, MRCAs = T,
#'                 exposures = T, drivers=T, add_driver_label = T)
my_annotate_forest <- function(tree_plot, forest, samples = TRUE, MRCAs = TRUE,
                            exposures = FALSE, facet_signatures = TRUE,
                            color_map = FALSE,
                            drivers = TRUE, add_driver_label = TRUE) {
  
  samples_info <- forest$get_samples_info()
  # Sampling times
  if (samples) {
    
    max_Y <- max(tree_plot$data$y, na.rm = TRUE)
    
    tree_plot <- tree_plot +
      ggplot2::geom_hline(
        data = samples_info,
        aes(
          yintercept = max_Y - time
          # color = name
        ),
        color = "indianred3",
        linetype = "dashed",
        linewidth = .5
      )
  }
  
  # MRCAs
  if (MRCAs) {
    sample_names <- samples_info %>% dplyr::pull(.data$name)
    
    MRCAs_cells <- lapply(sample_names,
                          function(s) {
                            forest$get_coalescent_cells(
                              forest$get_nodes() %>%
                                dplyr::filter(sample %in% s) %>%
                                dplyr::pull(.data$cell_id)
                            ) %>%
                              dplyr::mutate(sample = s)
                          }) %>%
      Reduce(f = dplyr::bind_rows) %>%
      dplyr::group_by(.data$cell_id) %>%
      dplyr::mutate(
        cell_id = paste(.data$cell_id)
      ) %>%
      dplyr::summarise(
        label = paste0("    ", .data$sample, collapse = "\n")
      )
    
    layout <- tree_plot$data %>%
      dplyr::select(.data$x, .data$y, .data$name) %>%
      dplyr::mutate(cell_id = paste(.data$name)) %>%
      dplyr::filter(.data$name %in% MRCAs_cells$cell_id) %>%
      dplyr::left_join(MRCAs_cells, by = "cell_id")
    
    
    tree_plot <-
      tree_plot +
      ggplot2::geom_point(
        data = layout,
        ggplot2::aes(x = .data$x, y = .data$y),
        color = "purple3",
        size = 3,
        pch = 21
      ) +
      ggplot2::geom_text(
        data = layout,
        ggplot2::aes(x = .data$x, y = .data$y,
                     label = .data$label),
        color = "purple3",
        size = 3,
        hjust = 0,
        vjust = 1
      )
  }
  
  if (exposures) {
    if (inherits(forest, "Rcpp_PhylogeneticForest")) {
      max_Y <- max(tree_plot$data$y, na.rm = TRUE)
      # Get exposures table
      exposures <- forest$get_exposures()
      
      exposure_colors <- get_colors_for(exposures %>%
                                          dplyr::pull(signature) %>%
                                          unique)
      
      # Add exposures start and end times for each signature
      times <- exposures$time %>%  unique() %>%  sort()
      
      exposures <- exposures %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
          t_end = dplyr::case_when(time == max(times) ~ Inf,
                                   .default = min(times[times >= time]))
        ) %>%
        dplyr::mutate(signature = factor(signature,
                                         levels = exposures %>%
                                           dplyr::arrange(time) %>%
                                           dplyr::pull(signature) %>%
                                           unique()))
      
      breaks <- sort(unique(exposures$exposure))
      
      # Annotate exposures on tree
      tree_plot <- tree_plot +
        ggplot2::geom_rect(
          data = exposures,
          ggplot2::aes(
            xmin = -Inf,
            xmax = Inf,
            ymin = ifelse(is.infinite(t_end), 0, max_Y - t_end),
            ymax = max_Y - time,
            fill = signature,
            alpha = exposure
          )
        ) +
        ggplot2::scale_fill_manual(values = exposure_colors) +
        ggplot2::scale_alpha_continuous(range = c(0.25, 0.75),
                                        breaks = breaks) +
        ggplot2::guides(fill = ggplot2::guide_legend(title = "Signature"),
                        alpha = ggplot2::guide_legend(title = "Exposure"))
      
      if (facet_signatures) {
        tree_plot <- tree_plot + ggplot2::facet_wrap(~ signature)
      }
      # Push exposure rectangles to the back
      layers_new <- list(tree_plot$layers[[length(tree_plot$layers)]])
      layers_new <- c(layers_new,
                      tree_plot$layers[1:(length(tree_plot$layers) - 1)])
      
      tree_plot$layers <- layers_new
    }
  }
  
  if (drivers) {
    if (inherits(forest, "Rcpp_PhylogeneticForest")) {
      drivers_mutations = drivers_CNAs = data.frame()
      try(expr = {
        drivers_codes <- forest$get_driver_mutations() %>% 
          dplyr::rename(from=start) %>% 
          dplyr::mutate(type=case_when(type=="SID"~"SNV", TRUE~type))
        drivers_mutations <- forest$get_sampled_cell_mutations() %>%
          dplyr::filter(class == "driver") %>%
          dplyr::left_join(drivers_codes,by=c("chr","from","ref","alt","type")) %>% 
          dplyr::mutate(driver_type = type) %>% 
          dplyr::rename(driver_id = code) %>% 
          # dplyr::mutate(driver_id = paste0(chr, ":", chr_pos, ":",
          #                                  ref, ">", alt),
          #               driver_type = type) %>%
          dplyr::select(cell_id, driver_id, driver_type)
      })
      try(expr = {
        drivers_CNAs <- forest$get_sampled_cell_CNAs() %>%
          dplyr::filter(class == "driver") %>%
          dplyr::mutate(driver_id = paste0(chr, ":", begin, "-",
                                           end, ":", allele),
                        driver_type = "CNA") %>%
          dplyr::select(cell_id, driver_id, driver_type)
      })
      
      drivers <- dplyr::bind_rows(drivers_mutations, drivers_CNAs)
      
      drivers_start_nodes <- lapply(unique(drivers$driver_id), function(d) {
        nodes_with_driver <- drivers %>%
          dplyr::filter(driver_id == d) %>%
          dplyr::pull(cell_id)
        d_type <- drivers %>%
          dplyr::filter(driver_id == d) %>%
          dplyr::pull(driver_type) %>%
          unique()
        
        forest$get_coalescent_cells(nodes_with_driver) %>%
          dplyr::mutate(driver_id = d, driver_type = d_type)
      }) %>%
        dplyr::bind_rows() %>%
        dplyr::mutate(cell_id = as.character(cell_id)) %>%
        dplyr::group_by(cell_id) %>%
        dplyr::summarise(driver_id = paste0(driver_id, collapse = "\n"))
      
      layout <- tree_plot$data %>%
        dplyr::select(x, y, name,mutant) %>%
        dplyr::mutate(cell_id = paste(name), has_driver = TRUE) %>%
        dplyr::filter(name %in% drivers_start_nodes$cell_id) %>%
        dplyr::left_join(drivers_start_nodes, by = "cell_id")
      
      tree_plot <- tree_plot +
        ggplot2::geom_point(
          data = layout,
          ggplot2::aes(x = .data$x, y = .data$y),
          fill = "#d68910",
          color = "#FF000000",
          size = 2,
          pch = 21
        )
      
      if (add_driver_label) {
        nudge_x <- (max(tree_plot$data$x) - min(tree_plot$data$x)) * .15
        tree_plot <- tree_plot +
          ggrepel::geom_label_repel(
            data = layout,
            ggplot2::aes(x = .data$x, y = .data$y,
                         label = .data$driver_id,fill=.data$mutant),
            color = "black",
            size = 2.5,
            hjust = 0,
            nudge_x = -nudge_x,
            direction = "x",
            alpha=0.5
          )+
          scale_fill_manual(values=color_map)
      }
    }
  }
  
  tree_plot
}

samples_table <- function(snapshot=NULL, sample_forest,use_snapshot=F,sim=NULL) {
  if (use_snapshot){
    sim <- ProCESS::recover_simulation(snapshot)
  }
  
  #sample_forest <- load_samples_forest(forest)
  info = sim$get_samples_info() ## requested from either the simulation recovery or as saved table
  
  nodes = sample_forest$get_nodes()
  clones = nodes %>% 
    dplyr::filter(!is.na(sample)) %>% 
    dplyr::group_by(sample, mutant) %>% 
    dplyr::pull(mutant) %>% 
    unique()
  clones_of_origin = nodes %>%
    dplyr::filter(!is.na(sample)) %>% 
    dplyr::group_by(sample, mutant) %>% 
    # dplyr::mutate(mutant = gsub(" ", "_", mutant)) %>% 
    dplyr::count(mutant) %>% 
    tidyr::pivot_wider(values_from = n, names_from = mutant, values_fill = 0) %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(Sample_Type = sum(c_across(clones) == 0)) %>% 
    dplyr::mutate(Sample_Type = ifelse(Sample_Type == (length(clones)-1), "Monoclonal", "Polyclonal")) %>% 
    dplyr::mutate(Total_Cells = sum(c_across(clones), na.rm = T)) %>% 
    dplyr::mutate(across(all_of(clones), ~ round(.x/Total_Cells,2), .names = "{.col} proportion"))
  
  samples_tb = dplyr::full_join(info, clones_of_origin, by = join_by("name" == "sample")) %>% 
    dplyr::select(!c("xmin","xmax","ymin","ymax","id","tumour_cells","tumour_cells_in_bbox")) %>% 
    dplyr::rename(Sample_ID=name) %>% 
    dplyr::rename(Samping_Time=time) %>% 
    dplyr::mutate(Samping_Time=round(Samping_Time,2)) %>% 
    dplyr::arrange(Sample_ID)
  
  
  
  return(samples_tb)
}
