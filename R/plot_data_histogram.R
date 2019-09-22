#' Plot the data histogram for a patient.
#' 
#' @description Plot the CCF histogram of each one of the samples available for
#' a patient, facetting the sample variable.
#'
#' @param x A REVOLVER cohort object.
#' @param patient The id of a patient.
#' @param cex Cex of the plot.
#' @param ... Extra parameters, not used.
#'
#' @return A \code{ggplot} plot.
#' 
#' @export
#' 
#'
#' @examples
#' TODO
plot_data_histogram = function(x,
                               patient,
                               cex = 1,
                               ...)
{
  # Samples and CCF are all we use to make this plot
  samples = Samples(x, patient)
  values = Data(x, patient) %>%
    select(!!samples, id, cluster, is.clonal)
  
  # Melt everything and make the first plot
  CCF_values = values %>%
    select(!!samples, id, cluster) %>%
    reshape2::melt(id = c('id', 'cluster'))
  
  # Clusters
  Cluster_values = values %>%
    select(cluster, id) %>%
    reshape2::melt(id = 'id')

  ncluster = length(unique(Cluster_values$value))

   ggplot(CCF_values %>%
            filter(value > 0), aes(value)) +
     geom_histogram(binwidth = 0.01) +
     facet_wrap(~variable) +
     theme_minimal(base_size = 10 * cex) +
     theme(
       legend.position = 'bottom',
       legend.key.size = unit(3, 'mm')
     ) +
     labs(
       title = "Data histogram",
       x = 'CCF',
       y = 'Observations'
     )
}
