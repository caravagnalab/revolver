#' Plot the dendrogram of REVOLVER"s clusters.
#' 
#' @description 
#' 
#' Plot the dendrogram of REVOLVER"s clusters, where leaves are patient
#' ids, coloured by cluster.
#'
#' @param x A \code{REVOLVER} object with fits and clusters.
#' @param clusters_palette A palette function that should return the colour of
#' an arbitrary number of clusters.
#' 
#' @family Plotting functions
#' 
#' @return A \code{ggplot} plot.
#' 
#' @import ggdendro
#' @export
#'
#' @examples
#' # Data released in the 'evoverse.datasets'
#' data('TRACERx_NEJM_2017_REVOLVER', package = 'evoverse.datasets')
#' 
#' plot_dendrogram(TRACERx_NEJM_2017_REVOLVER)
plot_dendrogram = function(x,
                           cluster_palette = distinct_palette_few)
{
  obj_has_clusters(x)
  
  # Dendrogram - it gives the ordering of the patients which are displayed on the x-axcis
  hc = x$cluster$fits$hc
  
  # Patient ordering - from the dedrogram, thi defines the levels of the factors used
  factors_patient_level = hc$order.lab
  
  # Get colors for the clusters
  clusters_colors = get_cluster_colors(x, cluster_palette)
  
  # Assign the colors following factors_patient_level
  patients_factors_colors = sapply(factors_patient_level,
                                   function(y)
                                     clusters_colors[Cluster(x, y) %>% pull(cluster)])
  names(patients_factors_colors) = factors_patient_level
  
  bars_separation = Cluster(x, factors_patient_level)
  bars_separation$cluster = factor(bars_separation$cluster, levels = unique(bars_separation$cluster))
  bars_separation = bars_separation %>% pull(cluster) %>% table %>% cumsum + 0.5
  
  # Number of clusters
  nclusters = Cluster(x) %>% pull(cluster) %>% unique %>% length
  
  # Dendrogram plot
  ggdendro::ggdendrogram(hc,
                         rotate = FALSE, 
                         size = 2) +
    my_ggplot_theme() +
    theme(axis.text.x = element_text(
      angle = 90,
      size = 8,
      color = patients_factors_colors
    )) +
    geom_vline(
      xintercept = bars_separation,
      size = .3,
      color = 'darkred',
      linetype = 'dashed'
    ) +
    labs(
      x = 'Patient',
      y = "REVOLVER evolutionary distance",
      title = 'REVOLVER cluster dendrogram',
      subtitle = x$annotation,
      caption = paste0('k = ', nclusters, ' clusters, ',
                       'n = ', x$n$patients, ' patients.')
    ) +
    scale_color_manual(values = clusters_colors) 
    
}
