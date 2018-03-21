
#' @title REVOLVER constructor for an object of type "rev_phylo": a tree.
#'
#' @details
#'
#' REVOLVER constructor to create an object of type "rev_phylo", which represents a model
#' for a patient. This is the same object regarldess data are CCF or binary values.
#'
#' @param M Input can be either an adjacency matrix, a dataframe with columns "from" and 'to' that
#' lists the edges of the model, or a list of edges (vector) in the format "A~B" to denote edge A --> B.
#' @param patient Patient ID for this tree.
#' @param dataset Data associated to this patient.
#' @param samples Thde IDs for the samples of this patient (e.g., "R1", "R2", ...)
#' @param score A score for this model, that we seek to maximize.
#' @param annotation A string to annotate this tree.
#'
#' @return An object of class \code{"rev_phylo"}
#' @export
#'
#' @import crayon
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
  M,
  patient,
  dataset,
  samples,
  score,
  annotation = '')
{
  add.GL = TRUE

  # Input can be either a matrix, a patient or a vector
  if(!is.matrix(M)) {
    if(is.data.frame(M)) M = DataFrameToMatrix(M)
    if(is.vector(M)) M = edgesToMatrix(M)
  }

  adj_mat = M

  # Dataset with certain columns
  required.cols = c('Misc', 'variantID', 'cluster', 'is.driver', 'is.clonal')
  if(!all(required.cols %in% colnames(dataset))) stop('Dataset should have the following columns:', required.cols)

  # Get samples and the table of the clusters
  samples_ID = samples

  CCF = clusters.table(dataset, samples_ID)
  data = binarize(dataset, samples_ID)

  patient_ID = patient
  nodes_ID = rownames(M)

  numRegions = nrow(data)
  numNodes = nrow(adj_mat)

  # Compute sample attachment
  attachments = colSums(data)
  attachments = attachments[attachments == 1]
  attachments = sapply(names(attachments),
    function(x){
      attachments[x] = rownames(data)[data[, x, drop = FALSE] == 1]
      })

  if(add.GL)
  {
    # special case for empty graphs with no edges
    if(sum(M) == 0) r = rownames(M)
    else r = root(adj_mat)

    adj_mat = rbind(MatrixToDataFrame(adj_mat), data.frame(from = 'GL', to = r))
    adj_mat = DataFrameToMatrix(adj_mat)
  }

  # Create a small template object which we can use with phylo-specific functions
  obj <-
    structure(
      list(
        adj_mat = adj_mat,
        dataset = dataset,
        root = root(adj_mat)
      ),
      class = "rev_phylo",
      call = match.call()
    )

  # Compute the information transfer
  transfer = information.transfer(obj)

  # Create the final output object
  obj <-
    structure(
      list(
        adj_mat = adj_mat,
        score = score,
        root = root(adj_mat),
        patient_ID = patient_ID,
        germline = add.GL,
        numRegions = numRegions,
        numNodes = numNodes,
        samples_ID = samples_ID,
        dataset = dataset,
        CCF = CCF,
        binary.data = data,
        nodes_ID = nodes_ID,
        attachments = attachments,
        annotation = annotation,
        transfer = transfer,
        modelname = paste("Phylogeny for patient", patient_ID)
        # REVOLVER_VERSION = REVOLVER_VERSION
      ),
      class = "rev_phylo",
      call = match.call()
    )

  return(obj)
}


#' @title Summary of a \code{"rev_phylo"} object.
#' @details  Reports some summary for a \code{"rev_phylo"} object, this is a bit more than just using print.
#'
#' @param x An object of class \code{"rev_phylo"}.
#' @param ... Extra parameters
#' @param print.stat Print or not statistics to screen.
#'
#' @return none
#' @export
#' @import crayon
#'
#' @examples
#' data(CRC.cohort)
#' summary(CRC.cohort$phylogenies[['adenoma_3']][[1]])
summary.rev_phylo <- function(x, ..., print.stat = TRUE){
  print.rev_phylo(x, ...)

  cat(cyan('\nCCF clusters:\n'))
  print(x$CCF)

  cat(cyan('\nDrivers:\n'))
  print(x$dataset[x$dataset$is.driver, ])

  cat(cyan('\nInformation Transfer:\n'))
  print(x$transfer)

  stat = stats.rev_phylo(x)

  if(print.stat)  {
    cat(yellow('\n[print.stat = T] Empirical probabilities:\n'))
    print.noquote(stat$probs)

    cat(yellow('\n[print.stat = T] Suppes\' conditions (binary data):\n'))
    print.noquote(stat$Suppes)

    # cat(yellow('\n[print.stat = T] Crossing rule (CCF):\n'))
    # print.noquote(stat$CCF.crossing.rule)

    cat(yellow('\n[print.stat = T] Pigeonhole principle (CCF):\n'))
    print.noquote(stat$CCF.pigeonhole)
  }

  cat(cyan('\nViolations:\n'))
  cat('Suppes\' Temporal Priority       : ', green(nrow(stat$Suppes) - stat$violations['tp']), red(stat$violations['tp']), '\n')
  # cat('\t\tCCF Crossing Rule       : ', green(nrow(stat$CCF.crossing.rule) - stat$violations['cr']), red(stat$violations['cr']), '\n')
  cat('Suppes\' Probability Raising     : ', green(nrow(stat$Suppes) - stat$violations['pr']), red(stat$violations['pr']), '\n')
  cat('CCF Pigenhole Principle         : ', green(nrow(stat$CCF.pigeonhole) * ncol(stat$CCF.pigeonhole) - stat$violations['pp']), red(stat$violations['pp']), '\n')

  cat(cyan('\nGoodness-of-fit: '), stat$gofit)
}

#' @title Print a \code{"rev_phylo"} object.
#' @details Print a summary for a \code{"rev_phylo"} object, which includes MSDOS-like console layout for trees.
#'
#' @param x An object of class \code{"rev_phylo"}.
#' @param ... Extra parameters
#' @param digits Digits to use.
#'
#' @return none
#' @export
#' @import crayon
#'
#' @examples
#' data(CRC.cohort)
#' CRC.cohort$phylogenies[['adenoma_3']][[1]]
print.rev_phylo <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  stopifnot(inherits(x, "rev_phylo"))

  M = x$adj_mat

  printPretty = function(node, indent, last)
  {
    cat(indent)
    if (last) {
      cat("\\-")
      indent = paste(indent, " ", sep ='')
    }
    else {
      cat("|-")
      indent = paste(indent, "| ", sep ='')
    }
    cat(node)


    if(node %in% names(x$attachments)) cat(paste(' [', x$attachments[node], ']', sep = ''))
    if(node != 'GL' &&  x$CCF[as.character(node), 'is.driver']) cat(
      sprintf(' :: %s',
              paste(driver(x, as.character(node)), collapse = ', '),
                # paste(x$dataset[x$dataset$is.driver & x$dataset$cluster == node, "variantID"],
                    collapse = ', '))
    cat('\n')

    cl = children(M, node)

    for(c in cl)
      printPretty(c, indent, c == cl[length(cl)])
  }

  cat(bgCyan('REVOLVER'))
  cat(blue(paste(' -- ', x$modelname, sep ='')), '\n')

  cat(cyan('\n\tNodes    :'), ncol(M), ifelse(x$germline, 'with GL', 'without GL'), '\n')
  cat(cyan('\tEdges    :'), sum(M), '\n')
  cat(cyan('\tSamples  :'), paste(x$samples_ID, collapse = ', '),  ifelse(x$numRegions > 1, '(Multi-region)', '(Single-sample)'),  '\n')

  # cat('Annotation:', x$annotation, '\n')

  cat(cyan('\n\tPhylogeny:\t\t\t'), sprintf('%s%s', cyan('Annotation: '), italic(x$annotation)),'\n')
  printPretty(x$root, "\t\t   ", TRUE)

  cat(cyan('\n\tTransfer:\n'))
  apply(x$transfer$drivers, 1, function(w) cat('\t\t  ', w[1], '--->', w[2], '\n'))


  cat(cyan('\n\tScore    :'),  green(x$score), '\n')


}



#' @title Compute statistics for a \code{"rev_phylo"} tree.
#'
#' @details
#' Compute some statistics for a \code{"rev_phylo"} tree. Thee include violations of the
#' pigeonhole principle, Suppes' conditions etc. Stats are compute no matter what
#' data is used to create the trees, but one should only consider those relevant
#' to the actual type of data.
#'
#' @param x An object of class \code{"rev_phylo"}.
#'
#' @return A list with all computed statistics.
#' @export
#'
#' @examples
#' data(CRC.cohort)
#' stats.rev_phylo(CRC.cohort$phylogenies[['adenoma_3']][[1]])
stats.rev_phylo = function(x)
{
  probs = function(b, i, j){
    b = rbind(b, GL = 0)
    mi = sum(b[, i])/nrow(b)
    mj = sum(b[, j])/nrow(b)
    joi = (b[, i] %*% b[, j]) / nrow(b)
    return(c(mi = mi, mj = mj, joi = joi))
  }

  M = MatrixToDataFrame(x$adj_mat)

  if(x$germline) M = M[M$from != 'GL', , drop = FALSE]

  if(nrow(M) == 0){
    return(
      list(
        probs = NULL,
        violations = NULL,
        gofit = 1,
        Suppes = NULL,
        CCF.crossing.rule = NULL,
        CCF.pigeonhole = NULL
      )
    )

  }

  p = data.frame(matrix(nrow = 0, ncol = 3), stringsAsFactors = FALSE)
  colnames(p) = c('p(i)', 'p(j)', 'p(i,j)' )

    p = apply(M, 1,
        function(e){
          probs(x$binary.data, e[1], e[2])
        })
    p = t(p)
    rownames(p) = apply(M, 1,function(e) paste(e[1], '->', e[2]))
    colnames(p) = c('p(i)', 'p(j)', 'p(i,j)' )

    TP = apply(
      p,
      1,
      function(e){
        if(e[1] > e[2]) return('>')
        if(e[1] == e[2]) return('=')
        return('<')
      })

    PR = apply(
      p,
      1,
      function(e){
        if(e[1] * e[2] < e[3]) return('>')
        if(e[1] * e[2] == e[3]) return('=')
        return('<')
      })

    df.Suppes = cbind(TP, PR)
    colnames(df.Suppes) = c('Temporal Priority', 'Probability Raising')

  # CCF data for this patient
  CCF = x$CCF[, x$samples_ID, drop = FALSE]

  # Direction penalty -- it was embedded in PP so it's now commented
  # direction.penalty = apply(
  #   M,
  #   1,
  #   function(e){
  #   CCF[e[1], , drop = FALSE] >= CCF[e[2], , drop = FALSE]
  # })
  #
  # if(ncol(CCF) == 1) direction.penalty = matrix(direction.penalty, ncol = 1)
  # else   direction.penalty = t(direction.penalty)
  # rownames(direction.penalty) = rownames(p)
  # colnames(direction.penalty) = x$samples_ID

  ## Branching penalty is the Pigeonhole Principle
  Mmatrix = DataFrameToMatrix(M)

  branching.penalty = sapply(
    x$nodes_ID,
    function(n){
      # print(n)
      cl = children(Mmatrix, n)
      if(length(cl) == 0) return(NULL)
      CCF[n, , drop = FALSE] >= colSums(CCF[cl, , drop = FALSE])
    })
  branching.penalty = branching.penalty[!unlist(lapply(branching.penalty, is.null))]
  branchp = Reduce(rbind, branching.penalty)
  rownames(branchp) = names(branching.penalty)

  # Scores summary -- we removed the direction penalty
  tp.viol = sum(as.integer(df.Suppes[, 1] == '<'))
  pr.viol = sum(as.integer(df.Suppes[, 2] == '<'))
  # cr.viol = sum(as.integer(!direction.penalty))
  ph.viol = sum(as.integer(!branchp))

  # violations = c(tp = tp.viol, pr = pr.viol, cr = cr.viol, pp = ph.viol)
  violations = c(tp = tp.viol, pr = pr.viol, pp = ph.viol)

  # Some alternative Goodness-Of-Fit measures. In the end, we take the
  # exact definition of the penalty for phylogenetic trees.
  # gofit = 2 * (x$numNodes - 1) * x$numRegions + x$numNodes  * x$numRegions
  # gofit = 1 - sum(violations)/gofit
  # gofit = nrow(branchp) * x$numRegions
  # gofit = 1 - ph.viol/gofit
  capture.output({ gofit = node.penalty.for.branching(list(M), CCF) })


  return(
    list(
      probs = p,
      violations = violations,
      gofit = gofit,
      Suppes = df.Suppes,
      # CCF.crossing.rule = direction.penalty, # Removed
      CCF.pigeonhole = branchp
    )
  )

}




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
#' @import RColorBrewer
#' @import igraph
#'
#' @examples
#' data(CRC.cohort)
#' plot(CRC.cohort$phylogenies[['adenoma_3']][[1]])
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

  nonw.colors = colorRampPalette(
    brewer.pal(
      palette,
      n = brewer.pal.info[palette, 'maxcolors'])) (length(nodesWithDriver))
  names(nonw.colors) = nodesWithDriver

  colors[nodesWithDriver] = nonw.colors

  for(c in x$nodes_ID) if(length(driver(x, c)) == 0) colors[c] = 'gainsboro'

  colors = c(colors, GL = 'white')

  add.alpha <- function(col, alpha=1){
    if(missing(col))
      stop("Please provide a vector of colors.")
    apply(sapply(col, col2rgb)/255, 2,
          function(x)
            rgb(x[1], x[2], x[3], alpha=alpha))
  }
  colors = add.alpha(colors, alpha)

  g = graph_from_adjacency_matrix(x$adj_mat)
  edList = as_edgelist(g)

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

  nodList = V(g)$name
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
       layout = layout.reingold.tilford(g, root = x$root),
       vertex.size = 20 * graph.cex,
       vertex.color = colors[V(g)$name],
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
