#' Plot the patients' jackknife trajectory frequencies.
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
#' @import ggrepel
#' @family Plotting functions
#' @export
#'
#' @examples
#' # Data released in the 'evoverse.datasets'
#' data('TRACERx_NEJM_2017_REVOLVER', package = 'evoverse.datasets')
#'
#' plot_jackknife_coclustering(TRACERx_NEJM_2017_REVOLVER)
plot_jackknife_trajectories_stability = function(x,
                                                 annotate_probability = .90,
                                                 annotate_percentage = .2
                                                 )
{
  obj_has_clusters(x)
  obj_has_jackknife(x)
  
  trajectories = Jackknife_trajectories_stability(x)
  
  relevant = trajectories %>%
    filter(prob_resamp > annotate_probability, num_patients > annotate_percentage)
  
  ggplot(trajectories, aes(x = prob_resamp, y = num_patients)) +
    geom_point(size = 1) +
    labs(
      x = "Probability across resamples",
      y = "Average percentage of patients with trajectory",
      title = "Jackknifed trajectories stability",
      subtitle = paste(x$annotation),
      caption = paste0('k = ', nclusters, ' clusters, ',
                       'n = ', x$n$patients, ' patients, ',
                       x$jackknife$params$resamples, ' resamples with leave ', x$jackknife$params$leave.out, ' out')
    ) +
    my_ggplot_theme() +
    xlim(0,1) +
    ylim(0,1) +
    ggrepel::geom_label_repel(
      data = relevant,
      aes(label = sprintf('%s \u2192 %s', from, to)),
      color = 'gray',
      fill = 'darkred',
      size = 3,
      xlim = c(NA, .5)
    )
}
