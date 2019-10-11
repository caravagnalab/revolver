#' Plot the DET index.
#'
#' @description
#' Plot the index of Divergent Evolutionary Trajectories, for a set of drivers
#' using function \code{DET_index}. The plot is a barplot with colours reflecting
#' the number of distinct incoming edges in each driver (species), and the height
#' reflecting the actual DET index value.
#'
#' @param x A \code{REVOLVER} object with fits.
#' @param ... Parmeters passed to function \code{DET_index} in order to compute
#' the index value.
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
#' plot_DET_index(TRACERx_NEJM_2017_REVOLVER)
#' 
#' # Passing parameters to DET_index
#' plot_DET_index(TRACERx_NEJM_2017_REVOLVER, min.occurrences = 5)
plot_DET_index = function(x,  ...)
{
  index = DET_index(x, ...)
  cex = 1
  
  print(index)
  
  index$driver = factor(index$driver, levels = index$driver)
  
  ggplot(index,
         aes(x = driver, y = DET_index, fill = N)
         )+
    geom_bar(stat = 'identity') +
    coord_flip() +
    scale_fill_distiller(palette = 'Spectral', direction = 1) +
    labs(
      title = 'DET index',
      y = "index",
      x = 'Driver',
      subtitle = x$annotation
      ) +
    guides(
      fill = guide_colorbar("Upstream drivers", barwidth  = unit(3 * cex, 'cm'))
    ) +
    my_ggplot_theme()
}
