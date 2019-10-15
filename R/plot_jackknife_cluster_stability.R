#' Plot the patients' jackknife cluster staability
#' 
#' @description 
#' 
#'  This function plots a barplot of the patients' jackknife cluster staability.
#'  
#'  The input vactor to this function can be obtained via \code{\link{Jackknife_cluster_stability}};
#'  in the visualisation the clusters are ordered by staability.
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
#' plot_jackknife_cluster_stability(TRACERx_NEJM_2017_REVOLVER)
plot_jackknife_cluster_stability = function(x,
                                  cluster_palette = distinct_palette_few)
{
  obj_has_clusters(x)
  obj_has_jackknife(x)
  
  stopifnot(palette %in% rownames(RColorBrewer::brewer.pal.info))
  
  data_barplot = Jackknife_cluster_stability(x) %>% 
    enframe(name = 'cluster') %>%
    rename(stability = value) %>%
    arrange(desc(stability))
    
  # Prepare factors for the clusters ordering and colours
  factors_level = Cluster(x) %>% arrange(cluster) %>% pull(patientID)
  
  clusters_colors = get_cluster_colors(x, cluster_palette)
  
  data_barplot$cluster = factor(data_barplot$cluster, levels = data_barplot$cluster)
  
  ggplot(data_barplot,
         aes(x = cluster, y = stability, fill = cluster)) +
    scale_fill_manual(values = clusters_colors) +
    geom_bar(stat = 'identity') +
    my_ggplot_theme() +
    coord_flip() +
    ylim(0, 1) +
    labs(
     x = 'Cluster',
     y = 'Stability',
     title = "Jackknifed cluster stability",
     subtitle = paste(x$annotation),
     caption = paste0('k = ', nclusters, ' clusters, ',
                      'n = ', x$n$patients, ' patients, ',
                      x$jackknife$params$resamples, ' resamples with leave ', x$jackknife$params$leave.out, ' out')
    )
}
