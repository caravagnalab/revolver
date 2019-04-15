

#' Title
#'
#' @param x
#' @param patient
#' @param rank
#' @param node_palette
#'
#' @return
#' @export
#'
#' @import ggraph
#' @import ggrepel
#'
#' @examples
plot_trees = function(x, patient, rank = 1, ...)
{
  if (!has_patient_trees(x, patient))
    stop("There are no trees for ", patient, ", aborting.")

  ntrees = length(Phylo(x, patient))

  lapply(
    rank,
    function(r)
    {
      plot.rev_phylo(
        x = Phylo(x, patient)[[r]],
        ...) +
        labs(caption = paste0("Ranks ", r, ' out of ', ntrees))
    }
  )
}
