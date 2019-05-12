#' Plot the DET index.
#'
#' @description
#' Plot the index of Divergent Evolutionary Trajectories, for a set of drivers
#' using function \code{DET_index}. The plot is a barplot with colours reflecting
#' the number of distinct incoming edges in each driver (species), and the height
#' reflecting the actual DET index value.
#'
#' @param x A REVOLVER cohort with fits.
#' @param cex Cex of the plot.
#' @param ... Parmeters passed to function \code{DET_index} in order to compute
#' the index value.
#'
#' @return A \code{ggplot} object for the plot
#' 
#' @export
#'
#' @examples
#' data(Breast.fit)
#' plot_DET_index(Breast.fit, min.occurrences = 5)
plot_DET_index = function(x, cex = 1, ...)
{
  index = DET_index(x, ...)
  
  print(index)
  
  index$driver = factor(index$driver, levels = index$driver)
  
  ggplot(index,
         aes(x = driver, y = DET_index, fill = N)
         )+
    geom_bar(stat = 'identity') +
    coord_flip() +
    scale_fill_distiller(palette = 'Spectral', direction = 1) +
    theme_minimal(base_size = 10 * cex) +
    labs(
      y = "DET index",
      x = 'Driver',
      title = x$annotation
      ) +
    guides(
      fill = guide_colorbar("Number of distinct incoming drivers", barwidth  = unit(3 * cex, 'cm'))
    ) +
    theme(
      legend.position = 'bottom',
      legend.key.size = unit(3 * cex, 'mm')
    ) 
}
