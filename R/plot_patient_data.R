#' Plot the data for a patient.
#' 
#' @description This function creates a complex plot for the
#' data of a patient, assembling plots returned from thw following
#' functions: \code{\link{plot_data_clusters}}, 
#' \code{\link{plot_data_clone_size}}, \code{\link{plot_data_mutation_burden}}
#' and \code{plot_patient_oncoprint}.
#'
#' @param x A REVOLVER cohort object.
#' @param patient A patient id.
#' @param ... Extra parameters passed to all the used plotting functions.
#'
#' @return A figure assembled with \code{ggpubr}.
#' @export
#'
#' @examples
#' TODO
plot_patient_data = function(x, patient, ...)
{
  require(ggpubr)
  
  # first panel
  first = ggarrange(
    plot_CCF_clusters(Phylo(x, patient, rank = 1)),
    plot_clone_size(Phylo(x, patient, rank = 1)),
    plot_patient_mutation_burden(x, patient),
    nrow = 1,
    ncol = 3
  )

  # second panel
  second = ggarrange(
    plot_patient_oncoprint(x, patient, ...),
    plot_CCF_histogram(x, patient, ...),
    ncol = 2,
    nrow = 1,
    widths = c(1, .3)
  )
  
  # Assembly
  ggarrange(
    first,
    second,
    ncol = 1,
    nrow = 2
  )
}

