#' Plot the data for a patient.
#'
#' @description
#'
#' This function creates a complex plot for the
#' data of a patient, assembling plots returned from the following
#' functions: 1) \code{\link{plot_data_clusters}},
#' 2) \code{\link{plot_data_clone_size}}, 3) \code{\link{plot_data_mutation_burden}}
#' and 4) \code{plot_patient_oncoprint}.
#'
#' @param x A \code{REVOLVER} cohort.
#' @param patient A patient id.
#' @param ... Extra parameters passed to all the used plotting functions.
#'
#' @family Plotting functions
#'
#' @return A figure assembled with \code{ggpubr}.
#' @export
#'
#' @import ggpubr
#' @import ctree
#'
#' @examples
#' # Data released in the 'evoverse.datasets'
#' data('TRACERx_NEJM_2017_REVOLVER', package = 'evoverse.datasets')
#'
#' plot_patient_data(TRACERx_NEJM_2017_REVOLVER, 'CRUK0001')
#' plot_patient_data(TRACERx_NEJM_2017_REVOLVER, 'CRUK0002')
plot_patient_data = function(x, patient, ...)
{
  # first panel
  first = ggarrange(
    plot_patient_CCF_clusters(Phylo(x, patient, rank = 1)),
    plot_clone_size(Phylo(x, patient, rank = 1)),
    plot_patient_mutation_burden(x, patient),
    nrow = 1,
    ncol = 3
  )

  # second panel
  second = ggarrange(
    plot_patient_oncoprint(x, patient, ...),
    plot_patient_CCF_histogram(x, patient, ...),
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

