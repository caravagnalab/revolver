#' @title REVOLVER constructor for a tree
#'
#' @details
#'
#' This constructor creates an object of type "rev_phylo", which represents a model
#' for a patient. This is the same type of object regarldess data are CCF or binary values.
#'
#' @param x A REVOLVER cohort.
#' @param patient The id of this patient in the cohort.
#' @param M The tree adjacency matrix.
#' @param score A score for this model, that we seek to maximize.
#' @param annotation A string to annotate this tree.
#'
#' @return An object of class \code{"rev_phylo"} that represents this tree.
#' @export
#'
#' @import crayon
#' @import tidygraph
#'
#' @examples
#' # Take data from CRC, inspect sample and dataset
#' data("CRC")
#' dataset = CRC[CRC$patientID == 'adenoma_3', ]
#' samples = paste('R', 1:5, sep = '')
#'
#' # Extract CCF values
#' CCF = sapply(dataset$CCF, revolver:::CCF.parser)
#' CCF = t(apply(CCF, 2, as.numeric))
#' rownames(CCF) = rownames(dataset)
#' colnames(CCF) = samples
#'
#' # Bind a dataset with explicit data
#' dataset = cbind(dataset, CCF)
#'
#' # Create empty adj_mat
#' m = matrix(0, ncol = 4, nrow = 4)
#' colnames(m) = rownames(m) = 1:4

#' # Phylogeny: will attach a GL (germline) to the root
#' # of the tree. In this case all nodes are roots.
#' revolver_phylogeny(m, 'adenoma_1', dataset, samples, 52)
revolver_phylogeny = function(
  x,
  patient,
  M,
  score,
  annotation = paste0("A tree for patient ", patient))
{
  # This function will create this output object
  obj <-
    structure(
      list(
        adj_mat = NA,        # Adjacency matrix for the tree
        tb_adj_mat = NA,     # Tidygraph object for this tree
        score = NA,          # Tree score
        patient = NA,        # Patient ID
        samples = NA,        # Samples name
        drivers = NA,        # Driver mapping to clusters
        CCF = NA,            # CCF (aggregated per cluster)
        annotation = NA,     # Some custom annotation
        transfer = NA,       # Information transfer from this tree
        tree_type = NA       # Mutation tree or phylogeny
      ),
      class = "rev_phylo",
      call = match.call()
    )

  # =-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-
  # The information that we need for each tree are
  # CCF clusters, data and driver events (mapped to clusters)
  # =-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-
  obj$CCF = CCF_clusters(x, patient)
  obj$samples = Samples(x, patient)
  obj$drivers = Drivers(x, patient) %>%
    select(variantID, cluster)

  # =-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-
  # We begin to create a representation of the tree
  # =-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-

  # Input can should be an adjacency matrix (which works also for monoclonal tumours)
  if(class(M) != 'matrix') stop("Input `M` should be an adjacency matrix, aborting.")

  adj_mat =  M
  df_mat = MatrixToDataFrame(M)

  # We chack for that to be a tree - empty is OK in this casem if monoclonal
  if(!is_tree(adj_mat, allow.empty = nrow(obj$CCF) == 1)) stop("The input adjacency matrix is not a valid tree, aborting.")

  # Add GL node, beware of the special case of empty adj_mat (might happen for a monoclonal tumour)
  M_root = ifelse(sum(adj_mat) == 0, rownames(adj_mat), root(adj_mat))
  df_mat = rbind(df_mat, data.frame(from = 'GL', to = M_root, stringsAsFactors = FALSE))

  # Update the adj matrix with GL, which can now go into obj
  obj$adj_mat = DataFrameToMatrix(df_mat)

  # =-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-
  # Create a tidygraph object for this tree
  # =-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-
  tb_adj_mat = as_tbl_graph(df_mat) %>%
    activate(nodes) %>%
    rename(cluster = name)

  # We can add information specific for this tree to the tidygraph
  tb_adj_mat = tb_adj_mat %>%
    left_join(obj$CCF, by = 'cluster')

  # Sample attachment for input data
  attachment = obj$CCF %>%
    select(cluster, obj$samples) %>%
    reshape2::melt(id = 'cluster') %>%
    mutate(value = ifelse(value > 0, 1, 0)) %>%
    group_by(variable) %>%
    filter(sum(value) == 1) %>%
    ungroup() %>%
    filter(value == 1) %>%
    select(-value) %>%
    rename(attachment = variable) %>%
    mutate(attachment = paste(attachment)) %>%
    group_by(cluster) %>%
    summarise(attachment = paste(attachment, collapse = ', '))

  tb_adj_mat = tb_adj_mat %>%
    left_join(attachment, by = 'cluster')

  # Drivers per node
  tb_adj_mat = tb_adj_mat %>%
    left_join(
      obj$drivers %>%
        group_by(cluster) %>%
        summarise(driver = paste(variantID, collapse = ', ')),
      by = 'cluster')

  # Store it in obj
  obj$tb_adj_mat = tb_adj_mat

  # Compute the information transfer
  obj$transfer = information_transfer(obj)

  # Extras
  obj$score = score
  obj$patient = patient

  obj$annotation = annotation

  obj$tree_type = ifelse(
    all(obj$CCF %>% select(obj$samples) %in% c(0, 1)),
    "Mutation tree (binary data)",
    "Phylogenetic tree (CCF)"
  )

  return(obj)
}


#' @title Print a \code{"rev_phylo"} object.
#' @details Print a summary for a \code{"rev_phylo"} object, which includes MSDOS-like console layout for trees.
#'
#' @param x An object of class \code{"rev_phylo"}.
#' @param ... Extra parameters
#'
#' @return Nothing
#'
#' @export
#'
#' @import crayon
#'
#' @examples
#' data(CRC.cohort)
#' CRC.cohort$phylogenies[['adenoma_3']][[1]]
print.rev_phylo <- function(x, ...)
{
  stopifnot(inherits(x, "rev_phylo"))

  M = x$adj_mat
  tb = x$tb_adj_mat

  printPretty = function(node, indent, last)
  {
    cat(indent)
    if (last) {
      cat("\\-")
      indent = paste(indent, " ", sep = '')
    }
    else {
      cat("|-")
      indent = paste(indent, "| ", sep = '')
    }
    cat(node)

    A = tb %>%
      activate(nodes) %>%
      filter(cluster == !!node) %>%
      pull(attachment)

    if (!is.na(A))
      cat(paste(' [', A, ']', sep = ''))

    D = tb %>%
      activate(nodes) %>%
      filter(cluster == !!node) %>%
      pull(driver)

    if (!is.na(D))
      cat(sprintf(' :: %s', D))

    cat('\n')

    cl = children(M, node)

    for (c in cl)
      printPretty(c, indent, c == cl[length(cl)])
  }

  pio::pioHdr(paste0('REVOLVER tree - ', x$annotation), suffix = '\n')

  cat('\n')
  print(x$CCF)

  pio::pioStr("Tree shape (drivers annotated)", "", prefix = '\n', suffix = '\n\n')

  printPretty(node = "GL", indent = "  ", last = TRUE)

  pio::pioStr("Information transfer", "", prefix = '\n', suffix = '\n\n')

  apply(x$transfer$drivers, 1, function(w)
    cat('  ', w[1], '--->', w[2], '\n'))

  pio::pioStr("Tree score", x$score, suffix = '\n', prefix = '\n')
}



#' S3 method that plots a REVOLVER tree.
#'
#' @param x A REVOLVER tree (object of class \code{"rev_phylo"}).
#' @param cex Cex for the plot.
#' @param node_palette A function that applied to a number will return a set of colors.
#' By default this is a \code{colorRampPalette} applied to 9 colours of the \code{RColorBrewer}
#' palette \code{Set1}. Colors are generated following a topological sort of the information
#' transfer, which is obtained from \code{igraph}.
#' @param tree_layout A layout that can be used by \code{tidygraph}, which wraps \code{igraph}'s
#' layouts. By default this is a `tree` layout.
#' @param add_information_transfer If `TRUE`, it will combine the plot of the input tree with
#' a plot for the associated information transfer for the driver events annotated. The colouring
#' of the nodes of the trees will match the colouring of the drivers. Combinations of plots is
#' done via \code{ggarrange} of package \code{ggpubr}.
#' @param icon If `TRUE` the icon tree version of a tree is plot. This type of view does not show 
#' a clone unless it has a driver annotated.
#' @param ... Extra parameters
#'
#' @return The plot. If `add_information_transfer = TRUE` the object is a combined figure from
#' package \code{ggpubr}, otherwise it is a single \code{ggplot} object.
#'
#' @export plot.rev_phylo
#'
#' @import crayon
#' @import igraph
#' @import tidygraph
#' @import ggraph
#' @import ggpubr
#' @import RColorBrewer
#'
#' @examples
#' data(CRC.cohort)
#' plot(CRC.cohort$phylogenies[['adenoma_3']][[1]])
plot.rev_phylo = function(x,
                          cex = 1,
                          node_palette = colorRampPalette(RColorBrewer::brewer.pal(n = 9, "Set1")),
                          tree_layout = 'tree',
                          add_information_transfer = FALSE,
                          icon = FALSE,
                          ...
                          )
{
  if(!icon)
  {
    # Get the tidygraph
    tree = x
    tb_tree = tree$tb_adj_mat
    
    # TODO Color edges as of information transfer
    #  - get path
    #  - modify edges etc.
    # tree$transfer
    
    # Color the nodes by cluster id, using a topological sort
    # to pick the colors in the order of appeareance in the tree
    clones_orderings = igraph::topo_sort(igraph::graph_from_adjacency_matrix(DataFrameToMatrix(tree$transfer$clones)),
                                         mode = 'out')$name
    
    nDrivers = length(clones_orderings) - 1 # avoid GL
    
    drivers_colors = c('white', node_palette(nDrivers))
    names(drivers_colors) = clones_orderings
    
    # Add non-driver nodes, with the same colour
    non_drivers = tb_tree %>%
      activate(nodes) %>%
      filter(!is.driver) %>%
      pull(cluster) # GL is not selected because is NA for is.driver
    
    non_drivers_colors = rep("gainsboro", length(non_drivers))
    names(non_drivers_colors) = non_drivers
    
    tb_node_colors = c(drivers_colors, non_drivers_colors)
    
    # Plot call
    layout <- create_layout(tb_tree, layout = tree_layout)
    
    mainplot = ggraph(layout) +
      geom_edge_link(
        arrow = arrow(length = unit(2 * cex, 'mm')),
        end_cap = circle(5 * cex, 'mm'),
        start_cap  = circle(5 * cex, 'mm')
      ) +
      geom_label_repel(
        aes(
          label = driver,
          x = x,
          y = y,
          colour = cluster
        ),
        na.rm = TRUE,
        nudge_x = .3,
        nudge_y = .3,
        size = 2.5 * cex
      ) +
      geom_node_point(aes(colour = cluster,
                          size = nMuts),
                      na.rm = TRUE) +
      geom_node_text(aes(label = cluster),
                     colour = 'black',
                     vjust = 0.4) +
      coord_cartesian(clip = 'off') +
      # theme_graph(base_size = 8 * cex, base_family = '') +
      theme_void(base_size = 8 * cex) +
      theme(legend.position = 'bottom',
            legend.key.size = unit(3 * cex, "mm")) +
      scale_color_manual(values = tb_node_colors) +
      scale_size(range = c(3, 10) * cex) +
      guides(color = FALSE,
             size = guide_legend(nrow = 1)) +
      labs(title = paste(tree$patient),
           subtitle = paste0('Scores ',
                             format(tree$score, scientific = T),
                             '.'))
    
    # Add information_transfer if required
    if (add_information_transfer)
    {
      mainplot = ggarrange(
        mainplot,
        plot_information_transfer(
          x,
          cex = cex,
          node_palette = node_palette,
          tree_layout = tree_layout,
          ...
        ),
        nrow = 1,
        ncol = 2
      )
    }
    
    return(mainplot)
  }
  else
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
}
