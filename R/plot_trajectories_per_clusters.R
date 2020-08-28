#' Plot the number of observed trajectories in each cluster.
#'
#' @description Assemble a figure that plots a graph (or tree) per
#' cluster where each driver is connected to its trajectories. This
#' allows to determine the most frequent trajectories in a cluster.
#'
#' @param x A revolver object with clusters computed.
#' @param min_counts A scalare >= 1 to subset only trajectories observed
#' in \code{min_counts} patients. If the value is in (0,1), the cut is
#' interpreted as a percentage and used to determine its actual value from
#' the cluster size (e.g., \code{n * min_counts}).
#'
#' @return A figure with multiple plots.
#' @export
#'
#' @examples
#' # Data released in the 'evoverse.datasets'
#' data('TRACERx_NEJM_2017_REVOLVER', package = 'evoverse.datasets')
#'
#' # In at least 5 cases
#' plot_trajectories_per_cluster(TRACERx_NEJM_2017_REVOLVER, min_counts = 5)
#'
#' # In 50% of the cluster's cases
#' plot_trajectories_per_cluster(TRACERx_NEJM_2017_REVOLVER, min_counts = .5)
plot_trajectories_per_cluster = function(x, min_counts = 5)
{
  stopifnot(revolver:::has_clusters(x))
  
  clusters_all = revolver::Cluster(x) %>% dplyr::pull(cluster) %>% unique() %>% sort()
  
  tree_p_plot = function(cl, min_counts = 1)
  {
    patients = revolver::Cluster(x) %>%
      dplyr::filter(cluster == cl) %>%
      dplyr::pull(patientID)
    
    revolver:::get_features(x, patients)
    
    # Get all trajectories
    All_trajectories = lapply(patients, function(p) {
      revolver::ITransfer(x,
                          p,
                          rank = 1,
                          type = "drivers",
                          data = "fits") %>%
        dplyr::mutate(patientID = p)
    })
    
    All_trajectories = Reduce(dplyr::bind_rows, All_trajectories) %>%
      dplyr::group_by(from, to) %>%
      dplyr::summarise(n = n()) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(desc(n)) %>%
      filter(n > min_counts)
    
    # Adapt min_counts
    if (min_counts > 0 & min_counts < 1) {
      cli::cli_alert_warning("min_counts in [0,1], interpreting that as a proportion.")
      min_counts = length(patients) * min_counts
    }
    
    All_trajectories = All_trajectories %>% filter(n > min_counts)
    
    if (nrow(All_trajectories) == 0)
      return(ggplot())
    
    G = tidygraph::as_tbl_graph(All_trajectories)
    
    ggraph(G, layout = 'tree') +
      geom_edge_diagonal(
        aes(
          start_cap = label_rect(from),
          end_cap = label_rect(to),
          edge_colour = n,
          # edge_width = n
        ),
        arrow = arrow(length = unit(2, 'mm')),
      ) +
      geom_node_text(aes(label = name),
                     colour = 'black',
                     vjust = 0.4) +
      coord_cartesian(clip = 'off') +
      theme_void(base_size = 10) +
      theme(legend.position = 'bottom',
            legend.key.width = unit(5, "mm")) +
      scale_edge_color_gradient(low = 'gray', high = 'indianred3') +
      # scale_edge_width(range = c(1, 3)) +
      guides(
        edge_color = guide_legend(title = "Occurrences", nrow = 1),
        edge_width = guide_legend(title = "Occurrences", nrow = 1)
      ) +
      labs(
        title = paste0('Cluster: ', cl, ' (n = ', length(patients), ")"),
        subtitle = paste0("Trajectories observed >", min_counts , " times"),
        caption = paste0(
          "Most/least frequent: ",
          All_trajectories$n[1],
          '/',
          All_trajectories$n[nrow(All_trajectories)]
        )
      )
    
  }
  
  ggpubr::ggarrange(plotlist = lapply(clusters_all, tree_p_plot, min_counts = min_counts))
  
}