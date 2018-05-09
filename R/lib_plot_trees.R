

#' Plot a phylogenetic or a mutation tree for a patient.
#'
#' @param x An object of class \code{"rev_phylo"}.
#' @param file Output file, or \code{NA}.
#' @param edge.width Edge width
#' @param edge.label Edge label
#' @param palette RColorBrewer palette to colour clusters.
#' @param graph.cex Cex for the graph.
#' @param table.cex Cex for the table.
#' @param alpha Transparency.
#' @param plot.stat \code{TRUE} if you want to show some statistics for the tree.
#' @param verbose Output.
#' @param ... Extra parameters
#'
#' @return Nothing
#' @export
#' @import crayon
#'
#' @examples
#' data(CRC.cohort)
#' plot(CRC.cohort$phylogenies[['adenoma_3']][[1]])
#'
#'


plot.rev_phylo = function(
  x, file = NA, edge.width = NA, edge.label = NA, palette = 'Set1',
  graph.cex = 1, table.cex = 1, alpha = 0.7, plot.stat = FALSE,
  verbose = FALSE, ...)
{
  # Get stats
  stat = stats.rev_phylo(x)
  if(is.null(stat$CCF.pigeonhole)) plot.stat = FALSE

  # Output devices etc.
  if(verbose) print(summary(x))
  if(!is.na(file)) pdf(file = file, ...)
  if(plot.stat) layout(matrix(c(1, 1, 2, 1,  1, 3), ncol = 2))

  colors = rep('white', x$numNodes)
  names(colors) = sort(x$nodes_ID)

  nodesWithDriver = sapply(x$nodes_ID, function(y) length(driver(x, y)) > 0)
  nodesWithDriver = names(nodesWithDriver)[nodesWithDriver]
  nodesWithDriver = sort(nodesWithDriver)

  nonw.colors = scols(nodesWithDriver, palette, alpha)
  names(nonw.colors) = nodesWithDriver

  colors[nodesWithDriver] = nonw.colors

  for(c in x$nodes_ID) if(length(driver(x, c)) == 0) colors[c] = 'gainsboro'

  colors = c(colors, GL = 'white')

  g = igraph::graph_from_adjacency_matrix(x$adj_mat)
  edList = igraph::as_edgelist(g)

  if(!all(is.na(edge.label))) {
    edLabel = apply(
      edList,
      1,
      function(x)
        if(x[1] == 'GL') paste('')
      else  paste(edge.label[x[1], x[2]]) )
  } else edLabel = ''


  if(!all(is.na(edge.width))) {
    edWeights = apply(
      edList,
      1,
      function(x) edge.width[x[1], x[2]] )
  } else edWeights = 1

  nodList = igraph::V(g)$name
  drivers = x$dataset[x$dataset$is.driver, ]

  drvCol = sapply(
    nodList,
    function(n) {
      if(n %in% drivers$cluster) return('red')
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

  for(i in 1:nrow(TR))
  {
    # Get the path
    p = find.path(x$adj_mat, from = TR[i,1], to = TR[i, 2])

    # color of the path
    path.color = colors[TR[i,1]]
    if(TR[i,1] == 'GL') path.color = 'darkblue'

    c = rep(path.color, nrow(p))
    names(c) = DataFrameToEdges(p)

    edColor.path = c(edColor.path, c)
  }

  for(i in 1:length(edColor.path))
  {
    edges = edColor.path[i]
    edColor[names(edges)] = edges
  }

  plot(g,
       layout = igraph::layout.reingold.tilford(g, root = x$root),
       vertex.size = 20 * graph.cex,
       vertex.color = colors[igraph::V(g)$name],
       # vertex.frame.color = drvCol,
       vertex.frame.color = 'white',
       edge.label = edLabel,
       edge.width =  edWeights * graph.cex * 2,
       edge.arrow.size = .5 * graph.cex,
       edge.color = edColor
  )

  title(bquote(bold('Patient: ') ~ .(x$patient_ID)),
        sub = x$annotation)

  legend('topright', title = 'Violations',
         legend = as.expression(
           c(
             # bquote(epsilon[cr] == ~ .(stat$violations['cr'])),
             bquote(epsilon[pp] == ~ .(stat$violations['pp'])),
             bquote(epsilon[tp] == ~ .(stat$violations['tp'])),
             bquote(epsilon[pr] == ~ .(stat$violations['pr']))
           )),
         col = 'darkred',
         pch = 19,
         bg = add.alpha('darkred', .2),
         box.col = 'white')

  legend('topleft',  title = 'Score',
         legend = as.expression(
           c(
             bquote(italic(f) == ~ .(x$score)),
             bquote(italic(g) == ~ .(stat$gofit))
           )),
         pch = 19,
         col = c('steelblue'),
         bg = add.alpha('steelblue', .3),
         box.col = 'white'
  )


  if(plot.stat)
  {
    colors = sapply(c('red', 'orange', 'darkgreen'), add.alpha, alpha = .8)
    names(colors) = c(-1, 0, 1)

    # Suppes' table
    stat$Suppes[stat$Suppes == '<'] = -1
    stat$Suppes[stat$Suppes == '='] = 0
    stat$Suppes[stat$Suppes == '>'] = 1

    table = matrix(apply(stat$Suppes, 2, as.numeric), ncol = nrow(stat$Suppes))
    colnames(table) = gsub(pattern = '->', '%->%', rownames(stat$Suppes))
    rownames(table) = c('TP', 'PR')

    image(table, col = colors[as.character(sort(unique(table)))], axes = FALSE)

    axis(2, at = seq(0, 1, length.out = ncol(table) ),
          labels = as.expression(parse(text = colnames(table))), las = 2 )
    axis(1, at = seq(0, 1, length.out = nrow(table) ),
         labels = rownames(table), las = 2 )

    title(main = expression(bold("Suppes\' conditions: ") * phantom("invalid (<), eq.(=), valid (>)")), col.main = 'black')
    title(main = expression(phantom(bold("Suppes\' conditions: ")) *
                              " invalid (<)," * phantom(", eq. (=), valid (>)")), col.main = colors[1])
    title(main = expression(phantom(bold("Suppes\' conditions: invalid (<),")) *
                              " eq. (=)," * phantom(", valid (>)")), col.main = colors[2])
    title(main = expression(phantom(bold("Suppes\' conditions: invalid (<), eq. (=),")) *
                              " valid (>)"), col.main = colors[3])

    # Pigeonhole principle
    table = matrix(apply(stat$CCF.pigeonhole, 2, as.numeric), ncol = nrow(stat$CCF.pigeonhole))
    rownames(table) = colnames(stat$CCF.pigeonhole)
    colnames(table) = rownames(stat$CCF.pigeonhole)

    image(table, col = colors[as.character(sort(unique(table)))], axes = FALSE)

    axis(2, at = seq(0, 1, length.out = ncol(table) ),
         labels = as.expression(parse(text = colnames(table))), las = 2 )
    axis(1, at = seq(0, 1, length.out = nrow(table) ),
         labels = rownames(table), las = 2 )

    title(main = expression(bold("Pigeonhole principle: ") * phantom("invalid (<=), valid (>)")), col.main = 'black')
    title(main = expression(phantom(bold("Pigeonhole principle: ")) *
                              " invalid (<=)," * phantom(", valid (>)")), col.main = colors[1])
    title(main = expression(phantom(bold("Pigeonhole principle: invalid (<=),")) *
                              " valid (>)"), col.main = colors[2])
  }

  if(!is.na(file)) dev.off()
}

# a-la-dictionary
driver = function(x, c){ x$dataset[which(x$dataset$cluster == c & x$dataset$is.driver), 'variantID'] }
driver.indexOf = function(x, d){ x$dataset[which(x$dataset$variantID == d & x$dataset$is.driver), 'cluster'] }

# Compute the frontier of "var" in a model. The frontier is the set of
# mutations in Gamma that are directly reachable via
# the transitive closure of the --> relation. So, these are the events
# selected by "var"'s evolutionary trajectories.
information.transfer = function(x, transitive.closure = FALSE, indistinguishable = FALSE)
{
  aux = function(r)
  {
    r.d = driver(x, r)

    if(length(r.d) > 0) return(r) # stop recursion

    c = children(model, r)

    if(is.null(c)) return(NULL) # leaf

    # recursion, reduction etc.
    return(
      Reduce(union,
             sapply(c, aux))
    )
  }

  model = x$adj_mat
  nodes.drivers = as.character(unique(x$dataset[x$dataset$is.driver, 'cluster']))

  # Root
  df = expand.grid(from = x$root, to = aux(x$root), stringsAsFactors = FALSE)

  # and all other stuff
  for(n in nodes.drivers)
  {
    df.n = NULL

    for(c in children(model, n))
      df.n = rbind(df.n, expand.grid(from = n, to = aux(c), stringsAsFactors = FALSE))

    df = rbind(df, df.n)
  }

  # Then, for all clones, we expand the drivers that they have
  expanded = apply(
    df,
    1,
    function(w) {
      e.from = driver(x, w['from'])
      if(length(e.from) == 0) e.from = w['from']

      e.to = driver(x, w['to'])

      expand.grid(from = e.from, to = e.to, stringsAsFactors = FALSE)
    })
  expanded = Reduce(rbind, expanded)

  # Then we see if we want to expand also the indistinguishible
  if(indistinguishable){
    for(n in nodes.drivers){
      expanded = rbind(expanded,
        expand.grid(from = driver(x, n), to = driver(x, n), stringsAsFactors = FALSE))
    }
  }

  # If we need transitive closures, we compute them here
  if(transitive.closure)
  {
    df = DataFrameToMatrix(df)
    df = nem::transitive.closure(df, mat = T, loops = FALSE)

    expanded = DataFrameToMatrix(expanded)
    expanded = nem::transitive.closure(expanded, mat = T, loops = FALSE)

    df = MatrixToDataFrame(df)
    expanded = MatrixToDataFrame(expanded)
  }

  return(list(clones = df, drivers = expanded))
}
