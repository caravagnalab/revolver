#' Plot REVOLVER trees for a patient
#'
#' @description
#'
#' This function is like \code{plot_data}, as it uses base plotting
#' functions to assemble a summary plot for a patient. This function
#' assembled the plot via \code{ggpubr}. The first line of plots
#' represents the top tree for a patient, and its information transfer.
#' The strip below represents up to the top-10 trees for this patient,
#' as they are obtained from the standard tree-scoring (which means
#' that the score is not accounting for the actual transfer, but just
#' for the tree structure).
#' 
#' @param x A REVOLVER cohort object
#' @param patient The patient for which the trees should be plot
#' 
#' @param ... Extra parameters, not used.
#'
#' @return A figure assembled with \code{ggpubr}.
#'
#' @export
#'
#' @import ggraph
#'
#' @examples
#' TODO
plot_trees = function(x, patient, ...)
{  
  require(ggpubr)
  
  if(!has_patient_trees(x, patient))
    stop(patient, " does not have the patient trees.")
  
  # =-=-=-=-=-=-=-=-
  # Top panel: top-ranking tree and its information transfer
  # =-=-=-=-=-=-=-=-
  
  top_tree = plot(Phylo(x, patient, rank = 1))
  top_tree_IT = plot(Phylo(x, patient, rank = 1), information_transfer = TRUE)
  
  first = ggarrange(
    top_tree,
    top_tree_IT,
    ncol = 2, 
    nrow = 1
  )
  
  # =-=-=-=-=-=-=-=-
  # Second panel: top-10 ranking trees in icon format
  # =-=-=-=-=-=-=-=-
  
  # Other 10 trees in 1x10 icon format
  ntrees =  min(
    Stats_trees(x, patient)$numTrees,
    10)
  
  # Strip plot
  strip = lapply(1:ntrees, Phylo, x = x, p = patient)
  
  strip = ggarrange(
    plotlist = lapply(strip, plot, icon = TRUE), 
    nrow = 1, ncol = 10)
  
  second = annotate_figure(
    strip,
    fig.lab = paste0("Top-",  ntrees, ' trees'))
  
  # Assembly
  ggarrange(
    first,
    second,
    heights = c(1, .3),
    ncol = 1,
    nrow = 2
  )
}

# plot_trees = function(x, patient, rank = 1, data = 'trees', ...)
# {
#   if (!has_patient_trees(x, patient))
#     stop("There are no trees for ", patient, ", aborting.")
# 
#   ntrees = length(Phylo(x, patient))
# 
#   lapply(
#     rank,
#     function(r)
#     {
#       tree_pl = plot.rev_phylo(x = Phylo(x, patient, r, data = data), ...) 
#       return(tree_pl)
#       
#       # if(data == 'trees')
#       #   return(
#       #     tree_pl +
#       #       labs(caption = paste0("Ranks ", r, ' out of ', ntrees))
#       #   )
#       # 
#       # if(data == 'fits')
#       #   return(
#       #     tree_pl +
#       #       labs(caption = paste0("Tree fit."))
#       #   )
#     }
#   )
# }
