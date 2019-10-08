#' Title
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
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
