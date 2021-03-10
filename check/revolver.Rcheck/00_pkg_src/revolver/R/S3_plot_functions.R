

#' Summary scatterplot of a cohort.
#'
#' @details Returns a scatterplot of the tumour mutational burden,  at the
#' clonal and subclonal level. Each dot is sized by the number of
#' clones with drivers, and coloured by the numnber of drivers.
#'
#' @param x A REVOLVER cohort.
#' @param ... Standard S3 signature.
#'
#' @family S3 functions
#'
#' @return A \code{ggplot} object of the plot.
#' @export
#' @examples
#' # Data released in the 'evoverse.datasets'
#' data('TRACERx_NEJM_2017_REVOLVER', package = 'evoverse.datasets')
#'
#' # Cancel the fits - just set the field to NULL
#' TRACERx_NEJM_2017_REVOLVER$fit = NULL
#'
#' plot(TRACERx_NEJM_2017_REVOLVER)
plot.rev_cohort = function(x,  ...)
{
  cex = 1

  ggplot(
    Stats(x),
    aes(
      color = numDriverMutations,
      size = numClonesWithDriver,
      x = numTruncalMutations,
      y = numSubclonalMutations
    )
  ) +
    stat_density_2d(size = .2 * cex, color = 'gray') +
    geom_point(alpha = .8) +
    geom_rug(size = .3 * cex) +
    scale_color_distiller(palette = 'Spectral') +
    scale_y_log10() +
    scale_x_log10() +
    coord_cartesian(clip = 'off') +
    labs(
      title = x$annotation,
      subtitle = paste("Cohort summary"),
      x = "Truncal mutations",
      y = "Subclonal mutations"
    ) +
    guides(
      color = guide_colorbar("Drivers", barwidth  = unit(3 * cex, 'cm')),
      size = guide_legend("Clones with drivers", barwidth  = unit(3 * cex, 'cm'))
    ) +
    my_ggplot_theme()

}


#' Summary scatterplot of a cohort's fits.
#'
#' @details Returns a scatterplot of the penalty of a patient's best fit,
#' versus the tumour mutational burden. Each dot is sized by the number of
#' combinations of Information Transfer for a patient, and coloured by the
#' number of annotated drivers.
#'
#' @param x A REVOLVER cohort with fits.
#' @param ... Standard S3 signature.
#'
#' @return A \code{ggplot} object of the plot.
#' @family S3 functions
#' @export 
#' @examples
#' # Data released in the 'evoverse.datasets'
#' data('TRACERx_NEJM_2017_REVOLVER', package = 'evoverse.datasets')
#'
#' plot(TRACERx_NEJM_2017_REVOLVER)
plot.rev_cohort_fit = function(x,  ...)
{
  obj_has_fit(x)
  cex = 1

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
    labs(
      subtitle = x$annotation,
      title = paste("Fits summary"),
      x = "(1 - penalty)",
      y = "Mutational burden"
    ) +
    guides(
      color = guide_colorbar("Drivers", barwidth  = unit(3 * cex, 'cm')),
      size = guide_legend("Transfers", barwidth  = unit(3 * cex, 'cm'))
    ) +
    scale_y_log10() +
    scale_x_continuous(limits = c(-0.1, 1.1)) +
    my_ggplot_theme()


  # Alternative plot - maybe for the future
  # icons = lapply(x$fit$phylogenies, plot_icon)
  # k = ceiling(sqrt(length(icons)))
  #
  # icons_figure = ggpubr::ggarrange(
  #   plotlist = icons,
  #   nrow = k,
  #   ncol = k
  # )
  #
  # icons_figure = ggpubr::annotate_figure(
  #   icons_figure,
  #   top = paste0("Trees (icon format) ~ ", x$annotation)
  # )
}
