#' Plot the occurrence of drivers across the cohort.
#' 
#' @description 
#' 
#' Plots a heatmap of the patients (x-axis) and driver events (y-axis),
#' which helps visualizing in a simple fashion the mapping of drivers
#' in the cohort.
#'
#' @param x A \code{REVOLVER} cohort.
#'
#' @return A `ggplot` object of the plot.
#' 
#' @family Plotting functions
#' 
#' @export
#'
#' @examples
#' # Data released in the 'evoverse.datasets'
#' data('TRACERx_NEJM_2017_REVOLVER', package = 'evoverse.datasets')
#'  
#' plot_drivers_occurrence(TRACERx_NEJM_2017_REVOLVER)
plot_drivers_occurrence = function(x)
{
  occurrences = lapply(
    x$patients, 
    Drivers, x = x)
  
  occurrences = Reduce(bind_rows, occurrences)
  
  occurrences = occurrences %>%
    select(patientID, variantID) 
  
  ggplot(
    occurrences,
    aes(
      y = variantID,
      x = patientID
    )
  ) +
    geom_tile() +
    my_ggplot_theme() +
    labs(      subtitle = x$annotation,
               title = "Drivers occurrence",
               x = 'Patient',
               y = 'Driver') +
    theme(
      axis.text.x = element_text(angle = 90),
    ) 
    
  
}
