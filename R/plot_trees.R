#' Plot REVOLVER trees for a patient
#'
#' @description
#'
#' This function plots one or more trees for a patient, in two
#' possible formats (icon and non-icon). The icon format is a tiny
#' representation of a tree with coloured nodes and edges, which
#' shows only the information transfer associated to the tree. The
#' non-icon version is the canonical representaion of the full tree,
#' with nodes coloured by drivers status, node size scaled by the clone
#' size (per node), and drivers annotated as labels. Switching from
#' one representation to the other is possible via a parameter
#' \code{icon}. This function returns a list of trees (which can
#' be assembled via \code{ggpubr} functions for instance); the number
#' of trees that one wants to plot is controlled via parameter \code{rank}.
#'
#'
#' @param x A REVOLVER cohort object
#' @param patient The patient for which the trees should be plot
#' @param rank A vector that specifies which tree should be plot for a
#' patient. The trees are stored according to the rank, and this is an 
#' integer vector which will be used to access the trees. By default, only
#' the top-rank tree is plot (\code{rank = 1}).
#' 
#' @param data 
#' @param ... Extra parameters which will be passed to all the 
#' delegate plotting calls  other functions.
#'
#' @return A list of plots.
#'
#' @export
#'
#' @import ggraph
#' @import ggrepel
#'
#' @examples
plot_trees = function(x, patient, rank = 1, data = 'trees', ...)
{
  if (!has_patient_trees(x, patient))
    stop("There are no trees for ", patient, ", aborting.")

  ntrees = length(Phylo(x, patient))

  lapply(
    rank,
    function(r)
    {
      tree_pl = plot.rev_phylo(x = Phylo(x, patient, rank, data = data), ...) 
      
      if(data == 'trees')
        return(
          tree_pl +
            labs(caption = paste0("Ranks ", r, ' out of ', ntrees))
        )
      
      if(data == 'fits')
        return(
          tree_pl +
            labs(caption = paste0("Tree fit."))
        )
    }
  )
}
