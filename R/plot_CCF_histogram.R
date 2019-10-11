#' Plot the data histogram for a patient.
#' 
#' @description 
#' 
#' If this cohort is working with CCF data, this function plots
#' the CCF histogram of each one of the samples available for
#' a patient, facetting the sample variable. 
#' 
#' Otherwise, this plot does not make much sense.
#'
#' @param x A REVOLVER cohort object.
#' @param patient The id of a patient.
#'
#' @return A \code{ggplot} plot.
#' 
#' @family Plotting functions
#' 
#' @export
#' 
#' @examples
#' # Data released in the 'evoverse.datasets'
#' data('TRACERx_NEJM_2017_REVOLVER', package = 'evoverse.datasets')
#'  
#' plot_CCF_histogram(TRACERx_NEJM_2017_REVOLVER, 'CRUK0001')
#' 
#' plot_CCF_histogram(TRACERx_NEJM_2017_REVOLVER, 'CRUK0002')
plot_CCF_histogram = function(x, patient)
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
     my_ggplot_theme() +
     labs(
       title = "Data histogram",
       x = 'CCF',
       y = 'Observations'
     )
}
