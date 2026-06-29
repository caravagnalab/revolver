# Declare tidyverse column-name variables used in NSE (non-standard evaluation)
# to suppress R CMD CHECK "no visible binding" NOTEs.
utils::globalVariables(c(
  ".", "IT", "Misc", "N", "N_tot", "Solution", "Var1", "Var2",
  "clonal", "cluster", "cluster_size", "colorRampPalette",
  "combInfTransf", "data", "diversity", "driver", "driverID", "edges",
  "from", "id", "is.clonal", "is.driver", "label", "mean_CCF", "nMuts",
  "name", "newSolution", "node", "node1.driver", "node2.driver", "nodes",
  "numBiopsies", "numClonal", "numClonesWithDriver", "numDriverMutations",
  "numMutations", "numSubclonal", "numSubclonalMutations",
  "numTruncalMutations", "num_patients", "occurrences", "p.value",
  "p_clonal", "p_subclonal", "p_tot", "patientID", "penalty",
  "prob_resamp", "psign", "score", "stability", "sumCCF", "sumCCF_cl",
  "to", "trajectory", "type", "value", "variable", "variantID",
  "fisher.test", "median", "str"
))

#' @import ggplot2
#' @importFrom dplyr `%>%` arrange bind_rows count desc distinct do filter full_join group_by left_join mutate n pull rename row_number select summarise summarize ungroup
#' @importFrom tidyr nest separate spread
#' @importFrom tibble as_tibble enframe
#' @importFrom grDevices colorRampPalette
#' @importFrom grid arrow unit
#' @importFrom stats as.dendrogram as.dist as.hclust fisher.test median order.dendrogram sd
#' @importFrom utils capture.output combn data str
#' @importFrom ggpubr annotate_figure ggarrange text_grob
#' @importFrom tidygraph activate as_tbl_graph
#' @importFrom ctree plot_CCF_clusters plot_clone_size plot_icon plot_information_transfer
#' @import evoverse.datasets
NULL

.onLoad <- function(libname, pkgname) {
  options(pio.string_fg_colour = crayon::bgYellow$black)
  invisible()
}

.onAttach <- function(libname, pkgname) {
  revolver_welcome_message = getOption('revolver_welcome_message', default = TRUE)
  if (revolver_welcome_message) {
    packageStartupMessage(
      'Loading revolver, \'Repeated Evolution in Cancer\'. Support : https://caravagn.github.io/revolver/')
    options(revolver_welcome_message = FALSE)
  }
}
