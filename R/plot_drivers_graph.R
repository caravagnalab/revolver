#' Plot graph-alike summary statistics for the cohort drivers.
#'
#' @description
#' Plot a graph with driver genes and annotate with
#' different summary statistics for the trajectories that
#' involve the drivers. This visualisation shows the frequency
#' of the driver in the cohort (node size), the penalty for
#' each pair of odering (edge thickness), the significance for
#' the pair of orderings as of a Fisher test (edge coloring)
#' and the overall heterogeneity upstream a driver as of the
#' DET index (node coloring). This function has parameters
#' to subset the computation to a list of predefined drivers,
#' or drivers associated to trajectories with a minimum
#' recurrence in the fits.
#'
#' @param x A REVOLVER object with fits.
#' @param drivers The list of drivers to consider, all by default.
#' See also function \code{\link{plot_penalty}}.
#' @param min.occurrences The minimum number of occurrences for
#' a trajectory to be considered, zero by default. See also
#' function \code{\link{plot_penalty}}.
#' @param alpha_level The significance level for the enrichment Fisher test.
#' @param ... Extra parameters passed to the \code{create_layout} function
#' by \code{ggraph}. For instance, passing \code{algorithm = 'kk'} and
#' \code{layout = 'igraph'} the `igraph` layout `kk` will be adopted.
#'
#' @return A `ggplot` object of the plot.
#'
#' @family Plotting functions
#'
#' @import ggraph
#'
#' @export
#'
#' @examples
#' # Data released in the 'evoverse.datasets'
#' data('TRACERx_NEJM_2017_REVOLVER', package = 'evoverse.datasets')
#'
#' # Base plot, can be quite crowded
#' plot_drivers_graph(TRACERx_NEJM_2017_REVOLVER)
#'
#' # Reduce the number of nodes cutting off low-frequencies one
#' plot_drivers_graph(TRACERx_NEJM_2017_REVOLVER, min.occurrences = 5)
#'
#' # As above, but with a more stringent test
#' plot_drivers_graph(TRACERx_NEJM_2017_REVOLVER, min.occurrences = 5, alpha_level = 0.01)
plot_drivers_graph = function(x,
                              drivers = x$variantIDs.driver,
                              min.occurrences = 0,
                              alpha_level = 0.05,
                              ...
                              )
{
  lp = list(...)
  # print(lp)

  # Subset E to make computations, and create a graph
  E = x$fit$penalty %>%
    filter(to %in% drivers, count >= min.occurrences)

  driver_stats = Stats_drivers(x) %>%
    rename(driver = variantID)

  # Get tests for enrichment via Fisher
  tests = enrichment_test_incoming_edge(E, alpha_level = alpha_level)

  # Get the DET index
  index = DET_index(x,
                    drivers = drivers,
                    min.occurrences = min.occurrences)

  # Create a graph for the plot
  G = as_tbl_graph(E) %>%
    rename(driver = name) %>%
    left_join(driver_stats, by = 'driver')

  # We need the index of the nodes to re-index the tests results
  nodeIdx = G %>% activate(nodes) %>% as_tibble() %>% pull(driver)
  nodeIdx = (1:length(nodeIdx))
  names(nodeIdx) = G %>% as_tibble() %>% pull(driver)

  tests = tests %>%
    mutate(
      from_drv = from,
      to_drv = to,
      from = nodeIdx[from],
      to = nodeIdx[to]
    )

  G = G %>%
    activate(edges) %>%
    left_join(tests, by = c('from', 'to'))

  # The DET is easier because is a property of the nodes
  G = G %>%
    activate(nodes) %>%
    left_join(index, by = c('driver'))

  # Plot call
  if(length(lp) == 0)
    layout <- create_layout(G, algorithm = 'kk', layout = 'igraph')
  else
    layout <- create_layout(G, ...)

  ggraph(layout) +
    geom_edge_link(
      aes(
          start_cap = label_rect(node1.driver),
          end_cap = label_rect(node2.driver),
          edge_colour = psign,
          edge_width = penalty
        ),
      arrow = arrow(length = unit(2, 'mm')),
    ) +
    geom_node_point(
      aes(
        size = N_tot,
        color = DET_index
        )
    ) +
    geom_node_text(aes(label = driver),
                   colour = 'black',
                   vjust = 0.4) +
    coord_cartesian(clip = 'off') +
    my_ggplot_theme() +
    theme(
      legend.key.size = unit(2.5, "mm")
      ) +
    scale_edge_color_manual(values = c(`TRUE` = 'darkorange', `FALSE` = 'gray')) +
    scale_color_gradient(low = 'steelblue', high = 'darkred') +
    scale_size(range = c(3, 10)) +
    scale_edge_width(range = c(.1, 1.1)) +
    guides(
      color = guide_legend(title = "DET index", nrow = 1),
      edge_width = guide_legend(title = "Penalty", nrow = 1),
      edge_color = guide_legend(title = paste0("Significant enrichment at level ", alpha_level), nrow = 1),
      size = guide_legend(title = "Driver counts", nrow = 1)
      ) +
    labs(
      title = paste('Driver to driver trajectories')
    )
}



