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

  # We chack for that to be a tree
  if(!is_tree(adj_mat)) stop("The input adjacency matrix is not a valid tree, aborting.")

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
    mutate(attachment = paste(attachment, collapse = ', '))

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
#' @param x An object of class \code{"rev_phylo"}.
#' @param file Output file, or \code{NA}.
#' @param palette RColorBrewer palette to colour clusters.
#' @param cex Cex for the graph.
#' @param alpha Transparency.
#' @param verbose Output.
#' @param ... Extra parameters
#'
#' @return Nothing
#' @export plot.rev_phylo
#' @import crayon
#'
#' @examples
#' data(CRC.cohort)
#' plot(CRC.cohort$phylogenies[['adenoma_3']][[1]])
plot.rev_phylo = function(x,
                          cex = 1,
                          node_palette = colorRampPalette(RColorBrewer::brewer.pal(n = 9, "Set1"))
                          )
{
  # Get the tidygraph
  tree = x
  tb_tree = tree$tb_adj_mat

  # Color the nodes by cluster id
  nDrivers = sum(tb_tree %>% pull(is.driver), na.rm = TRUE)

  tb_node_colors = tb_tree %>% pull(cluster)
  names(tb_node_colors) = tb_node_colors

  tb_node_colors["GL"] = "white"
  tb_node_colors[!tb_tree %>% pull(is.driver)] = 'gainsboro'
  tb_node_colors[!is.na(tb_tree %>% pull(driver))] = node_palette(nDrivers)

  # Plot call
  layout <- create_layout(tb_tree, layout = 'dendrogram')

  ggraph(layout, 'dendrogram') +
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
    geom_node_point(
      aes(
        colour = cluster,
        size = nMuts
        ),
      na.rm = TRUE
    ) +
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
    labs(
      title = paste(tree$patient),
      subtitle = paste0(
        'Scores ',
        format(tree$score, scientific = T),
        '.'
      )
    )
}
