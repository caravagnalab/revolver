#' Plot a phylogenetic or a mutation tree for a patient.
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
#' @export
#' @import crayon
#'
#' @examples
#' data(CRC.cohort)
#' revolver_plt_tree(CRC.cohort$phylogenies[['adenoma_3']][[1]])
revolver_plt_tree = function(x,
                             type = 'CCF',
                             file = NA,
                             palette = 'Set1',
                             cex = 1,
                             alpha = 0.7)
{
  # Get stats
  stat = stats.rev_phylo(x)
  labels.stat = inlineStat(x, type)

  edge.width = NA
  edge.label = NA

  if (length(labels.stat) == 0)
    labels.stat = 'None'

  lot = mylayout.on(file, 1, c(4, 2), cex)

  colors = rep('white', x$numNodes)
  names(colors) = sort(x$nodes_ID)

  nodesWithDriver = sapply(x$nodes_ID, function(y)
    length(driver(x, y)) > 0)
  nodesWithDriver = names(nodesWithDriver)[nodesWithDriver]
  nodesWithDriver = sort(nodesWithDriver)

  nonw.colors = scols(nodesWithDriver, palette, alpha)
  names(nonw.colors) = nodesWithDriver

  colors[nodesWithDriver] = nonw.colors

  for (c in x$nodes_ID)
    if (length(driver(x, c)) == 0)
      colors[c] = 'gainsboro'

  colors = c(colors, GL = 'white')

  g = igraph::graph_from_adjacency_matrix(x$adj_mat)
  edList = igraph::as_edgelist(g)

  if (!all(is.na(edge.label))) {
    edLabel = apply(edList,
                    1,
                    function(x)
                      if (x[1] == 'GL')
                        paste('')
                    else
                      paste(edge.label[x[1], x[2]]))
  } else
    edLabel = ''


  if (!all(is.na(edge.width))) {
    edWeights = apply(edList,
                      1,
                      function(x)
                        edge.width[x[1], x[2]])
  } else
    edWeights = 1

  nodList = igraph::V(g)$name
  drivers = x$dataset[x$dataset$is.driver,]

  drvCol = sapply(nodList,
                  function(n) {
                    if (n %in% drivers$cluster)
                      return('red')
                    return('white')
                  })
  vertex.size = rep(1, x$numNodes)
  # print(drvCol)

  # Color paths based on information transfer
  edList.order = data.frame(from = edList[, 1], to = edList[, 2])
  edList.order = DataFrameToEdges(edList.order)

  edColor = rep('gainsboro', length(edList.order))
  names(edColor) = edList.order

  TR = x$transfer$clones
  edColor.path = NULL

  for (i in 1:nrow(TR))
  {
    # Get the path
    p = find.path(x$adj_mat, from = TR[i, 1], to = TR[i, 2])

    # color of the path
    path.color = colors[TR[i, 1]]
    if (TR[i, 1] == 'GL')
      path.color = 'darkblue'

    c = rep(path.color, nrow(p))
    names(c) = DataFrameToEdges(p)

    edColor.path = c(edColor.path, c)
  }

  for (i in 1:length(edColor.path))
  {
    edges = edColor.path[i]
    edColor[names(edges)] = edges
  }

  plot(
    g,
    layout = igraph::layout.reingold.tilford(g, root = x$root),
    vertex.size = 20,
    vertex.color = colors[igraph::V(g)$name],
    # vertex.frame.color = drvCol,
    vertex.frame.color = 'white',
    edge.label = edLabel,
    edge.width =  edWeights * 2,
    edge.arrow.size = .5,
    edge.color = edColor
  )

  title(bquote(bold('Patient: ') ~ .(x$patient_ID)),
        sub = x$annotation)


  legend(
    'topright',
    title = 'Violations',
    legend = labels.stat,
    col = 'darkred',
    pch = 19,
    bg = add.alpha('darkred', .2),
    box.col = 'white'
  )

  # legend(
  #   'topright',
  #   title = 'Violations',
  #   legend = as.expression(c(
  #     # bquote(epsilon[cr] == ~ .(stat$violations['cr'])),
  #     bquote(epsilon[pp] == ~ .(stat$violations['pp'])),
  #     bquote(epsilon[tp] == ~ .(stat$violations['tp'])),
  #     bquote(epsilon[pr] == ~ .(stat$violations['pr']))
  #   )),
  #   col = 'darkred',
  #   pch = 19,
  #   bg = add.alpha('darkred', .2),
  #   box.col = 'white'
  # )

  x$score = format(x$score, scientific = TRUE,  digits = 1)
  stat$gofit = format(stat$gofit, digits = 3)


  legend(
    'topleft',
    title = 'Score',
    legend = as.expression(c(
      bquote(italic(f) == ~ .(x$score)),
      bquote(italic(g) == ~ .(stat$gofit))
    )),
    pch = 19,
    col = c('steelblue'),
    bg = add.alpha('steelblue', .3),
    box.col = 'white'
  )


  mylayout.off(lot)
  invisible(NULL)
}




inlineStat = function(x, type = 'CCF')
{
  x = stats.rev_phylo(x)

  if (type == 'CCF')
    x = x$CCF.pigeonhole
  else
    x = x$Suppes

  if (is.null(x))
    return(list())

  dfedges = data.frame(stringsAsFactors = F)

  for (i in 1:nrow(x)) {
    for (j in 1:ncol(x)) {
      dfedges = rbind(
        dfedges,
        data.frame(
          entry.row = rownames(x)[i],
          entry.col = colnames(x)[j],
          value = x[i, j],
          stringsAsFactors = FALSE
        )
      )
    }
  }

  if (type == 'CCF')
    dfedges = dfedges[dfedges$value == FALSE, , drop = FALSE]
  else
    dfedges = dfedges[dfedges$value == '<', , drop = FALSE]


  string = apply(dfedges, 1, function(w)
    paste(w[1], w[2]))

  # return(paste(string, collapse = ', '))
  return(string)
}
