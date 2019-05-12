#' Print a REVOLVER cohort object with fits.
#'
#' @param x A REVOLVER cohort object from class \code{"rev_cohort_fit"}
#' @param ... unused arguments from generic
#'
#' @return Nothing
#' @export print.rev_cohort_fit
#'
#' @examples
#' data(Breast.fit)
#' Breast.fit
print.rev_cohort_fit = function(x, ...)
{
  class(x) = 'rev_cohort'
  print(x)
}


#' Fits summary scatterplot.
#' 
#' @details Returns a scatterplot of the penalty of a patient's best fit,
#' versus the tumour mutational burden. Each dot is sized by the number of 
#' combinations of Information Transfer for a patient, and coloured by the
#' number of annotated drivers. 
#'
#' @param x A REVOLVER cohort with fits.
#' @param cex Cex of the plot.
#'
#' @return A \code{ggplot2} object of the plot.
#' 
#' @export plot.rev_cohort_fit
#'
#' @examples
#' data(Breast.fit)
#' plot.rev_cohort_fit(Breast.fit)
plot.rev_cohort_fit = function(x, cex = 1, ...)
{
  obj_has_fit(x)
  
  TB = x$fit$fit_table
  
  TB = TB %>%
    left_join(Stats(x), by = "patientID")
  
  ggplot(TB,
         aes(
           x = penalty, 
           y = numMutations, 
           size = combInfTransf, 
           color = numDriverMutations
           )
         ) +
    stat_density_2d(size = .2 * cex, color = 'gray') +
    geom_point(alpha = .8) +
    geom_rug(size = .3 * cex) +
    scale_color_distiller(palette = 'Spectral') +
    theme_minimal(base_size = 10 * cex) +
    theme(
      legend.position = 'bottom',
      legend.key.size = unit(3 * cex, 'mm')
    ) +
    labs(
      title = x$annotation,
      subtitle = paste("Fits summary"),
      x = "(1 - penalty)",
      y = "Mutational burden"
    ) +
    guides(
      color = guide_colorbar("Drivers", barwidth  = unit(3 * cex, 'cm')),
      size = guide_legend("Transfers", barwidth  = unit(3 * cex, 'cm'))
    ) +
    scale_y_log10() +
    scale_x_continuous(limits = c(-0.1, 1.1)) 
}