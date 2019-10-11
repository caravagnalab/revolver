#' Plot the clonality status of cohort driver events.
#' 
#' @description 
#' 
#' Driver events can be annotated in clonal or subclonal clusters.
#' This function reports this information in a barplot.
#'
#' @param x A \code{REVOLVER} cohort.
#'
#' @return A \code{ggplot} object of the plot.
#' 
#' @family Plotting functions
#' 
#' @export
#'
#' @examples
#' # Data released in the 'evoverse.datasets'
#' data('TRACERx_NEJM_2017_REVOLVER', package = 'evoverse.datasets')
#'  
#' plot_drivers_clonality(TRACERx_NEJM_2017_REVOLVER)
plot_drivers_clonality = function(x)
{
  # Get drivers, mirrored
  st = Stats_drivers(x) %>%
    arrange(desc(numClonal), desc(numSubclonal))

  order_drivers = st$variantID
  st$numSubclonal = - st$numSubclonal

  st = st %>%
    select(variantID, numClonal, numSubclonal) %>%
    rename(Clonal = numClonal, Subclonal = numSubclonal) %>%
    reshape2::melt(id = 'variantID')

  st$variantID = factor(st$variantID, levels = order_drivers)

  # Percentages
  N = length(x$patients)

  ggplot(st,
         aes(x = variantID,
             y = value,
             fill = variable)) +
    geom_bar(stat="identity", position="identity") +
    coord_flip(clip = 'off') +
    scale_fill_manual(values = c(`Clonal` = 'steelblue', `Subclonal` = 'darkorange3')) +
    labs(
      title = paste("Driver burden"),
      y = paste0('Occurrences (n = ', N, ' patients)'),
      x = 'Driver',
      subtitle = paste(x$annotation)
    ) +
    guides(
      fill = guide_legend("")
    ) +
    my_ggplot_theme()
}
