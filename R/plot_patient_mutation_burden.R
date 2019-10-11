#' Plot the mutation burden for a patient.
#' 
#' @description This function creates a scatterplot
#' of a full cohort annotating clonal versus subclonal
#' mutational burden, colouring each point by the 
#' number of total drivers in the patient, and using
#' as shape the number of clones with drivers.
#'
#' @param x A REVOLVER cohort object.
#' @param patient A patient id.
#' @param cex Cex of the plot.
#' @param ... Extra parameters, not used.
#'
#' @family Plotting functions
#' 
#' @return A \code{ggplot} plot.
#' 
#' @export
#'
#' @examples
#' # Data released in the 'evoverse.datasets'
#' data('TRACERx_NEJM_2017_REVOLVER', package = 'evoverse.datasets')
#'  
#' plot_patient_mutation_burden(TRACERx_NEJM_2017_REVOLVER, 'CRUK0001')
#' plot_patient_mutation_burden(TRACERx_NEJM_2017_REVOLVER, 'CRUK0002')
plot_patient_mutation_burden = function(x, patient, ...)
{
  p = patient
  cex = 1
  st = Stats(x)
  
  # Clonal Subclonal
  ggplot(
    st,
    aes(
      x = numTruncalMutations,
      y = numSubclonalMutations,
      color = numDriverMutations,
      shape = factor(numClonesWithDriver)
    )
  ) +
    stat_density_2d(alpha = .3, 
                    aes(      
                      x = numTruncalMutations,
                      y = numSubclonalMutations
                      ),
                    color = 'black', inherit.aes = FALSE) +
    geom_point(size = 2 * cex) +
    ggrepel::geom_label_repel(
      data = st %>% filter(patientID == !!p),
      aes(label = patientID),
      nudge_x = .1,
      nudge_y = .1, 
      size = 3 * cex, 
      min.segment.length = 0
    ) +
    coord_cartesian(clip = 'off') +
    scale_color_gradient(low = 'steelblue', high = 'darkred') +
    labs(
      title = paste("Mutational burden for ", patient),
      x = "Clonal burden",
      y = "Subclonal burden",
      subtitle = paste0(
        st %>% filter(patientID == !!patient) %>% pull(numBiopsies), " sample(s) and ",
        st %>% filter(patientID == !!patient) %>% pull(numMutations), " mutation(s)."
        )
    ) +
    scale_x_log10() +
    scale_y_log10() +
    guides(
      color = guide_colourbar("Drivers", barwidth = 6),
      shape = guide_legend("Clones with drivers")
    ) +
    my_ggplot_theme()
    
}