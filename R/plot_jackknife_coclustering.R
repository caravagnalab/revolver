#' Plot the patients' jackknife co-clustering probability
#' 
#' @description 
#' 
#'  This function plots a raster of the patients' jackknife co-clustering probability.
#'  
#'  The input matrix to this function can be obtained via \code{\link{Jackknife_patient_coclustering}};
#'  in the visualisation the patient ids are coloured and ordered according to the clusters in the
#'  input dataset, which can be obtained via \code{\link{Cluster}}.
#'
#' @param x A \code{REVOLVER} cohort with fits, clusters and jackknife results available.
#' @param cluster_palette A palette of colours; must be a function that when applied to
#' a number it returns that number of colours.
#'
#' @return A \code{ggplot} figure.
#' @family Plotting functions
#' @export
#'
#' @examples
#' # Data released in the 'evoverse.datasets'
#' data('TRACERx_NEJM_2017_REVOLVER', package = 'evoverse.datasets')
#'
#' plot_jackknife_coclustering(TRACERx_NEJM_2017_REVOLVER)
plot_jackknife_coclustering = function(x,
                                       cluster_palette = distinct_palette_few)
{
  obj_has_clusters(x)
  obj_has_jackknife(x)

  co_clustering = Jackknife_patient_coclustering(x)

  # Prepare factors for the clusters ordering and colours
  factors_level = Cluster(x) %>% arrange(cluster) %>% pull(patientID)

  nclusters = Cluster(x) %>% pull(cluster) %>% unique %>% length

  clusters_colors = get_cluster_colors(x, cluster_palette)
  
  factors_colors =  Cluster(x) %>%
    arrange(cluster) %>%
    mutate(color = clusters_colors[cluster])

  # Long data format and ordering
  longData = reshape2::melt(co_clustering) %>%
    filter(value > 0)

  longData$Var1 = factor(longData$Var1, levels = factors_level)
  longData$Var2 = factor(longData$Var2, levels = factors_level)

  ggplot(longData, aes(x = Var2, y = Var1)) +
    geom_raster(aes(fill = value)) +
    scale_fill_viridis_c(direction = -1) +
    labs(
      x = "Patients",
      y = "Patients",
      title = "Jackknifed patient co-clustering probability",
      subtitle = paste(x$annotation),
      caption = paste0('k = ', nclusters, ' clusters, ',
                       'n = ', x$n$patients, ' patients, ',
                       x$jackknife$params$resamples, ' resamples with leave ', x$jackknife$params$leave.out, ' out')
    ) +
    my_ggplot_theme() +
    guides(fill = guide_colourbar("Probability", barwidth = unit(3, 'cm'))) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1,
                                 color = factors_colors$color),
      axis.text.y = element_text(color = factors_colors$color)
    ) +
    geom_vline(
      xintercept = factors_colors %>% group_by(cluster) %>% summarise(n = n()) %>% pull(n) %>% cumsum + .5,
      size = .3,
      color = 'gray',
      linetype = 'dashed'
    ) +
    geom_hline(
      yintercept = factors_colors %>% group_by(cluster) %>% summarise(n = n()) %>% pull(n) %>% cumsum + .5,
      size = .3,
      color = 'gray',
      linetype = 'dashed'
    )
}
