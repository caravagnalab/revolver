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
#' @param rankA A vector that specifies which tree should be plot for a
#' patient. The trees are stored according to the rank. By default, only
#' the top-rank tree is plot (\code{rank = 1}).
#' @param node_palette A function that can return, for an input number,
#' a number of colours.
#'
#' @return A list of plots.
#'
#' @export
#'
#' @import ggraph
#' @import ggrepel
#'
#' @examples
plot_trees = function(x, patient, rank = 1, icon = FALSE, ...)
{
  if (!has_patient_trees(x, patient))
    stop("There are no trees for ", patient, ", aborting.")

  ntrees = length(Phylo(x, patient))

  lapply(
    rank,
    function(r)
    {
      # Mode 1: icon tree
      if(icon) {
        return(
          plot_tree_icon(x = Phylo(x, patient)[[r]], ...)
        )
      }

      # Mode 2: full tree
      plot.rev_phylo(
        x = Phylo(x, patient)[[r]],
        ...) +
        labs(caption = paste0("Ranks ", r, ' out of ', ntrees))
    }
  )
}

# This is a function that plots a simplified version of a tree (icon format)
plot_tree_icon = function(
  x,
  cex = 1,
  node_palette = colorRampPalette(RColorBrewer::brewer.pal(n = 9, "Set1")),
  tree.layout = 'tree',
  ...
  )
{
  # Get the tidygraph
  tree = x
  tb_tree = tree$tb_adj_mat

  # Color the nodes by cluster id
  tb_node_colors = tb_tree %>% filter(is.driver) %>% pull(cluster)

  tb_node_colors = node_palette(length(tb_node_colors))
  tb_node_colors = c(tb_node_colors, `GL` = 'white')
  names(tb_node_colors) = c(tb_tree %>% filter(is.driver) %>% pull(cluster), 'GL')

  # Graph from transfer
  tb_icon = as_tbl_graph(tree$transfer$clones) %>%
    rename(cluster = name) %>%
    activate(edges) %>%
    mutate(
      cluster = .N()$cluster[from]
    )

  # Plot call
  layout <- create_layout(tb_icon, layout = tree.layout)

  ggraph(layout) +
    geom_edge_link(
      aes(colour = cluster)
      ) +
    geom_node_point(
      aes(colour = cluster),
      na.rm = TRUE,
      size = 3
    ) +
    coord_cartesian(clip = 'off') +
    theme_void(base_size = 8 * cex) +
    theme(legend.position = 'none') +
    scale_color_manual(values = tb_node_colors) +
    scale_edge_color_manual(values = tb_node_colors) +
    guides(color = FALSE,
           size = guide_legend(nrow = 1)
           )
  }
