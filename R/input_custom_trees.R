#' Add a custom set of trees for a patient.
#'
#' @description
#'
#' For a patient it is possible to add a set of precomputed trees, which need
#' to be stored as adjacency matrices, and must be scored and ranked by scoring.
#' The trees must have as nodes all the clusters that are annotated in the data
#' of the patient. If the patient has already some trees, these are overwritten.
#'
#' @param x A REVOLVER cohort object.
#' @param patient A patient in the cohort, for which the trees are created.
#' @param trees A list a of precomputed adjacency matrices to describe the input trees.
#' @param scores A vectore of scores for the input trees.
#'
#' @return A modififed REVOLVER cohort with available phylogeneis for \code{patient}.
#'
#' @export
#' @family Cohort creation
#'
#' @import ctree
#'
#' @examples
#' # Data released in the 'evoverse.datasets'
#' data('TRACERx_NEJM_2017_REVOLVER', package = 'evoverse.datasets')
#'
#' # To make it simple we use some trees that are already available
#' trees = Phylo(TRACERx_NEJM_2017_REVOLVER, p = "CRUK0001")
#'
#' # Get matrices from these objects, and remove the GL columns/ rows
#' matrices = lapply(trees, function(x) x$adj_mat[rownames(x$adj_mat) != 'GL', colnames(x$adj_mat) != 'GL'])
#'
#' print(matrices)
#'
#' # Get scores for these matrices - vector format
#' scores = sapply(trees, function(x) x$score)
#'
#' print(scores)
#'
#' # Add these trees - this function overwrites the original ones
#' TRACERx_NEJM_2017_REVOLVER = input_custom_trees(
#'   TRACERx_NEJM_2017_REVOLVER,
#'   patient = "CRUK0001",
#'   trees = matrices,
#'   scores = scores
#' )
input_custom_trees = function(
  x,
  patient,
  trees,
  scores
)
{
  stop_not_revolver_object(x)
  stop_invalid_patient(x, patient)

  pioTit("Using custom trees for ", patient)

  are_suitable_precomputed_trees(x, patient,
                                 precomputed.trees = trees,
                                 scores)

  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Run ctree constructor, preparing the return object
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  if(!has_patient_trees(x)) x$phylogenies = NULL

  x$phylogenies[[patient]] =
    easypar::run(
      FUN = function(w){
        ctree::ctree(
          CCF_clusters(x, patient),
          Drivers(x, patient),
          Samples(x, patient),
          patient = patient,
          M = trees[[w]],
          score = scores[[w]]
        )
      },
      PARAMS = lapply(seq_along(trees), list),
      parallel = FALSE,
      progress_bar = FALSE
    )

  # Just show how many combinations we have
  comb = combination_of_information_transfer(x, patient)
  pio::pioStr('\n Combinations of Information Transfer : ', comb, suffix = '\n')

  return(x)
}


# Check format of precomputed trees
are_suitable_precomputed_trees = function(
  x,
  patient,
  precomputed.trees = NULL,
  precomputed.scores = NULL
)
{
  # N clusters
  required_clusters = CCF_clusters(x, patient)$cluster

  # Check if they're all N x N
  cmt_size = sapply(precomputed.trees, ncol)
  rmt_size = sapply(precomputed.trees, nrow)

  OK_size = (cmt_size == rmt_size) && (cmt_size == length(required_clusters))

  if(!all(OK_size)) {
    message("Input trees must be NxN adjacency matrices because this patient has N = ",
            length(required_clusters), " clusters. Some are not.")

    stop("Cannot use these trees, aborting")
  }

  # Check if they have all N clusters
  OK_clusters = sapply(
    precomputed.trees,
    function(w) {
      all(required_clusters %in% colnames(w)) &&
        all(required_clusters %in% rownames(w))
    }
  )

  if(!all(OK_clusters)) {
    message("Input trees must contain entries for all clusters ",
            paste(required_clusters, collapse = ', '), " but some are not.")

    stop("Cannot use these trees, aborting")
  }

  # Check if they're actual trees
  OK_trees = sapply(precomputed.trees, ctree:::is_tree)

  if(!all(OK_trees)) {
    message("Input trees are not trees and might have multiple roots, disconnected components etc.")

    stop("Cannot use these trees, aborting")
  }

  # All scores are not NA
  OK_scores_na = !(sapply(precomputed.scores, is.na))

  if(!all(OK_scores_na)) {
    message("Some trees are scored as NA, while they should all be real values.")

    stop("Cannot use these trees, aborting")
  }

  # Sorted scores
  OK_scores_sort = !(is.unsorted(rev(precomputed.scores)))

  if(!all(OK_scores_sort)) {
    message("Trees should be passed in descreasing order of scores, but they are not.")

    stop("Cannot use these trees, aborting")
  }
}
