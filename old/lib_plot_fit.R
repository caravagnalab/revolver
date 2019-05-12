#' Plot a barplot of the penalty per patient.
#'
#' @param x A \code{"rev_cohort_fit"} object
#' @param file Output file, leave it to NA to plot to canvas.
#'
#' @return none
#' @export
#'
#' @importFrom graphics title
#'
#' @examples
#' data(Breast.fit)
#' revolver_plt_penalty_barplot(Breast.fit)
revolver_plt_penalty_barplot = function(x, file = NA, cex = 1)
{
  obj_has_fit(x)
  pio::pioHdr(
    'REVOLVER Plot: barplot of penalty across patients ',
    toPrint = c(`Output file` = file),
    prefix = '\t'
  )

  lot = mylayout.on(file, 1, c(4, 2), cex)

  # Trivial penalty barplots
  pen = x$fit$penalty
  pen = pen[order(pen, decreasing = T)]

  barplot(
    pen,
    log = 'x',
    col = scols(names(x$fit$penalty), 'YlOrRd'),
    pch = 19,
    xlab = ,
    border = NA,
    ylab = 'Patient',
    horiz = T,
    las = 2,
    cex.names = .5
  )

  title(bquote(bold('Penalty: ') ~ 'log' ~ 'p(' ~ T[i] ~ '|' ~ bold(w) ~
                 ')' ^ alpha),
        sub = 'Value (the lower, the worst)')


  mylayout.off(lot)
  invisible(NULL)
}

#' Plot a pheatmap of the penalty (patients per alterations) .
#'
#' @description
#' This function plots the multinomial model "w". It can plot its
#' value after the first step of the fitting algorithm, or the second.
#' It plots normalized or unnormalized values.
#'
#' @param x A \code{"rev_cohort_fit"} object
#' @param palette RColorBrewer palette
#' @param normalized Set to TRUE to normalize the penalty.
#' @param type Use this to select either 'after.expansion', or 'before.expansion'
#' estimates of the penalty
#' @param file Output file, leave it to NA to plot to canvas.
#'
#' @return none
#' @export
#'
#' @examples
#' data(Breast.fit)
#' revolver_plt_penalty_matrix(Breast.fit)
revolver_plt_penalty_matrix = function(x,
                                       palette = 'Blues',
                                       normalized = FALSE,
                                       type = 'after.expansion',
                                       file = NA,
                                       cex = 1)
{
  obj_has_fit(x)
  pio::pioHdr(
    'REVOLVER Plot: matrix of penalty across patients and alterations',
    toPrint = c(
      `Normalize w` = normalized,
      `Use fits before/ after expansion` = type,
      `Output file` = file
    ),
    prefix = '\t'
  )

  stopifnot(type %in% c('after.expansion', 'before.expansion'))
  stopifnot(palette %in% rownames(RColorBrewer::brewer.pal.info))

  # Counts data
  counts = NULL

  if (type == 'before.expansion')
  {
    counts = x$fit$multinomial.penalty
  }
  else
  {
    feat = revolver.featureMatrix(x)
    counts = DFW2Matrix(feat$consensus.explosion)
  }

  # Ordering is the same
  counts = counts[order(colnames(counts)), order(colnames(counts))]

  # Colours and breaks
  colors = scols(1:max(counts), palette)
  breaks = -1:max(counts)

  # Update these if normalization is required
  if (normalized)
  {
    counts = sweep(counts, 2, colSums(counts), FUN = "/")
    colors = scols(seq(0.0001, 1.1, by = 0.01), palette)
    breaks = c(-0.00001, seq(0.0001, 1.1, 0.01))
  }

  pheatmap::pheatmap(
    counts,
    main = paste(
      "Penalty across patients/ alterations\n Normalized",
      normalized,
      'Type',
      type
    ),
    breaks = breaks,
    color = c('gainsboro', colors),
    fontsize_row = 6,
    border = NA,
    fontsize_col = 6,
    cellwidth = 8,
    cellheight = 8,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    file = file
  )

  invisible(NULL)
}


#' Plot the DET index, via non parametric bootstrap.
#'
#' @description
#' Plot the index of Divergent Evolutionary Trajectories, computed via
#' a non parametric bootstrap procedure. The current DET is implemented
#' in two different ways: via Shannon's entropy or the Variance Analog index.
#' You can also decide if you want to use models computed after  likelihood fit, or
#' after likelihood fit and expansion.
#'
#' @param x A \code{"rev_cohort_fit"} object
#' @param N Number of non-parametric bootstrap resamples.
#' @param DET.type Index type, default is \code{"Shannon"} (Shannon's equitability,
#' which is entropy), also available is \code{"VA"} (Variance Analog index in \code{VA})
#' @param type Use this to select either 'after.expansion', or 'before.expansion'
#' estimates of the penalty
#' @param colour Colour of the bars.
#' @param file Output file, or NA.
#'
#' @return The DET values
#' @export
#'
#' @importFrom graphics abline barplot hist
#'
#'
#' @examples
#' data(Breast.fit)
#' revolver_plt_DET_index(Breast.fit)
revolver_plt_DET_index = function(x,
                                  DET.type = 'Shannon',
                                  type = 'after.expansion',
                                  N = 100,
                                  colour = 'steelblue',
                                  file = NA,
                                  cex = 1)
{
  obj_has_fit(x)

  pio::pioHdr(
    'REVOLVER Plot: boostrapped DET index (full cohort)',
    toPrint = c(
      `Type of index` = DET.type,
      `Use fits before/ after expansion` = type,
      `Number of non parametric bootstrap replicates` = N,
      `Output file` = file
    ),
    prefix = '\t'
  )

  stopifnot(type %in% c('after.expansion', 'before.expansion'))
  stopifnot(DET.type %in% c('Shannon', 'VA'))

  lot = mylayout.on(file, 1, c(4, 2), cex)

  DET = revolver_DETindex(
    x,
    n.boot = N,
    table = T,
    index = DET.type,
    type = type
  )$DET.cohort

  min = min(DET) - 2
  max = max(DET) + 2

  DET.type = ifelse(DET.type == 'Shannon',
                'Shannon\'s equitability',
                'Variance Analog Index')

  type = ifelse(type == 'after.expansion',
                'Likelihood-fit/ Expansion',
                'Likelihood-fit')
  hist(
    DET,
    col = colour,
    breaks = 22,
    border = NA,
    main = bquote(bold('DET(w) full cohort') ~ 'N =' ~ .(N)),
    sub = bquote(.(DET.type) ~  italic(.(type))),
    ylab = NA,
    xlab = NA,
    xlim = c(min, max)
  )
  abline(
    v = median(DET),
    col = 'red',
    lty = 2,
    lwd = 2
  )

  mylayout.off(lot)
  return(DET)
}




#' Plot the DET index of each driver from the data (empirical).
#'
#' @description
#' The same index of \code{\link{revolver_plt_DET_index}} but computed per driver
#' and without the bootstrap.
#'
#' @param x A \code{"rev_cohort_fit"} object
#' @param DET.type Index type, default is \code{"Shannon"} (Shannon's equitability,
#' which is entropy), also available is \code{"VA"} (Variance Analog index in \code{VA})
#' @param type Use this to select either 'after.expansion', or 'before.expansion'
#' estimates of the penalty
#' @param palette RColorBrewer palette
#' @param file Output file, or NA.
#'
#' @return The DET values
#' @export
#'
#'
#' @examples
#' data(Breast.fit)
#' revolver_plt_DET_index_driver(Breast.fit)
revolver_plt_DET_index_driver = function(x,
                                         DET.type = 'Shannon',
                                         type = 'after.expansion',
                                         palette = 'Blues',
                                         file = NA,
                                         cex = 1)
{
  obj_has_fit(x)
  pio::pioHdr(
    'REVOLVER Plot: empirical DET index (per driver)',
    toPrint = c(
      `Type of index` = DET.type,
      `Use fits before/ after expansion` = type,
      `Output file` = file
    ),
    prefix = '\t'
  )

  stopifnot(type %in% c('after.expansion', 'before.expansion'))
  stopifnot(DET.type %in% c('Shannon', 'VA'))

  lot = mylayout.on(file, 1, c(4, 2), cex)

  DET = revolver_DETindex(
    x,
    n.boot = 1,
    table = T,
    index = DET.type,
    type = type,
    driver.only = TRUE
  )$DET.driver

  DET.type = ifelse(DET.type == 'Shannon',
                    'Shannon\'s equitability',
                    'Variance Analog Index')

    type = ifelse(type == 'after.expansion',
                'Likelihood-fit/ Expansion',
                'Likelihood-fit')

  barplot(
    sort(DET),
    col = scols(seq(0, 1, 0.01), palette),
    border = NA,
    pch = 19,
    xlab = NA,
    ylab = NA,
    main = bquote(bold('DET(w) per driver') ~ italic('empirical')),
    sub = bquote(.(DET.type) ~ '   ' ~ italic(.(type))),
    horiz = T,
    las = 2,
    cex.names = .5
  )

  mylayout.off(lot)
  return(DET)
}



#' Plot the tree fit for a patient.
#'
#' @description
#' This function uses REVOLVER's standard plotting functions but shows the fit
#' for a model.
#'
#' @param x A \code{"rev_cohort_fit"} object
#' @param patient The patient to plot.
#' @param palette RColorBrewer palette.
#' @param alpha Alpha colouring.
#' @param file Output file, or NA.
#' @param cex Cex for graphics.
#'
#' @return None
#' @export
#'
#'
#' @examples
#' data(Breast.fit)
#' revolver_plt_fit_patient(Breast.fit, patient = 'PD9850')
revolver_plt_fit_patient = function(x,
                                    patient,
                                    palette = 'Set1',
                                    alpha = .8,
                                    file = NA,
                                    cex = 1)
{
  obj_has_fit(x)

  pio::pioHdr(
    'REVOLVER Plot: fit for a patient',
    toPrint = NULL, prefix = '\t'
  )

  phylo = x$fit$phylogenies[[patient]]

  revolver_plt_tree(phylo, file = file, palette = palette, alpha = alpha, cex = cex)
}

#' @title Plot the trajectories in a patient's fit model.
#'
#' @description
#' From the phylogenetic tree fit to this patient, REVOLVER extracts and plots the
#' trajectories depicted by the tree.
#'
#' @param x A \code{"rev_cohort_fit"} object
#' @param patient The patient to plot.
#' @param palette RColorBrewer palette.
#' @param alpha Alpha colouring.
#' @param file Output file, or NA.
#' @param cex Cex for graphics.
#'
#' @return None
#' @export
#'
#'
#' @examples
#' data(Breast.fit)
#' revolver_plt_trajectories_patient(Breast.fit, patient = 'PD9850')
revolver_plt_trajectories_patient = function(x,
                                             patient,
                                             palette = 'Set1',
                                             alpha = .8,
                                             file = NA,
                                             cex = 1)
{
  obj_has_fit(x)

  # pio::pioHdr(
  #   paste('REVOLVER Plot: trajectories from the fit of', patient),
  #   toPrint = c(
  #     `Patient` = paste(patient, collapse = ', '),
  #     `Output file` = file
  #   ),
  #   prefix = '\t'
  # )

  pio::pioHdr(
    paste('REVOLVER Plot: trajectories from the fit of', patient),
    toPrint = NULL,
    prefix = '\t'
  )
  stopifnot(palette %in% rownames(RColorBrewer::brewer.pal.info))

  lot = mylayout.on(file, 1, c(4, 2), cex)

  ##### Extract the trajectories
  result = tryCatch({
    phylo = x$fit$phylogenies[[patient]]
    subst = x$fit$substitutions[[patient]]

    # get groups, and the adj matrix
    groups = unique(phylo$dataset[phylo$dataset$is.driver, 'cluster'])
    adj_matrix = x$fit$exploded[[patient]]

    G = igraph::graph_from_adjacency_matrix(adj_matrix)

    # visit groups in order according to a topological sort, ensures correct expansions
    TS = wrapTS(phylo$adj_mat)
    TS = TS[TS %in% groups]

    # colors for the nodes...
    igraph::V(G)$color  = "gainsboro"
    igraph::V(G)["GL"]$color <- "white"

    ncolors = length(subst)
    colors = scols(sort(names(subst)), palette)
    colors = add.alpha(colors, alpha)

    # nodes get colored according to their cluster color
    for (g in TS) {
      m = subst[[g]]

      for (n in colnames(m))
        igraph::V(G)[n]$color = colors[g]
    }

    # we work on the layout as well
    lay = igraph::layout.reingold.tilford(G, root = 'GL',  mode = 'all')
    rownames(lay) =  igraph::V(G)$name

    # This topological sort can work only if G has no loops
    if(igraph::is_dag(G))
    {
      TS = wrapTS(adj_matrix)
      for (node in TS)
        lay = fixLayer(adj_matrix, node, lay)
    }

    mark.groups = lapply(subst, function(w)
      colnames(w))
    mark.col = add.alpha(colors[names(subst)], alpha / 2)
    mark.border = colors[names(subst)]

    # print this plot, and merge it to the phylogeny
    plot(
      G,
      vertex.frame.color = 'white',
      edge.arrow.size = .5,
      mark.groups = mark.groups,
      mark.col = mark.col,
      mark.border = mark.border,
      edge.color = 'steelblue',
      # vertex.color = colors[membership(wc)],
      layout = lay
    )
    title(bquote(bold(.(patient) ~ ':') ~ "evolutionary trajectories"), sub = "Exploded trajectories")
  },
  warning = function(w) {
    cat(crayon::red('WARNING'), '| ')
    warning(w)
  },
  error = function(w) {
    cat(crayon::red('ERROR'), '| ')
    warning(w)
  },
  finally = {

  })


  mylayout.off(lot)
  invisible(NULL)
}



#' @title Plot the information transfer from a patient's fit.
#'
#' @details
#' From the phylogenetic tree fit to this patient, REVOLVER extracts and plots the
#' information transfer depicted by the tree. The tree is annotated with the number
#' of clonal and subclonal observations of each driver, and the amount of times each
#' edge is detected
#'
#' @param x A \code{"rev_cohort_fit"} object
#' @param patient The patients to plot.
#' @param file Output file, or NA.
#' @param cex Cex for graphics.
#'
#' @return None
#' @export
#' @import crayon
#'
#' @examples
#' data(Breast.fit)
#' revolver_plt_itransfer_patient(Breast.fit, patient = 'PD9850')
revolver_plt_itransfer_patient = function(x,
                                          patient,
                                          file = NA,
                                          cex = 1)
{
  obj_has_fit(x)

  pio::pioHdr(
    paste('REVOLVER Plot: information transfer from the fit of', patient),
    toPrint = NULL,
    prefix = '\t'
  )

  occ.table = clonal.subclonal.table(x)
  edge.table = rev_table.orderings(x, intelli.cutoff = 0)

  lot = mylayout.on(file, 1, c(4, 2), cex)

  ##### Extract the information transfer
  result = tryCatch({

    phylo = x$fit$phylogenies[[patient]]
    subst = x$fit$substitutions[[patient]]

    # get groups, and the adj matrix
    groups = unique(phylo$dataset[phylo$dataset$is.driver, 'cluster'])
    adj_matrix = x$fit$exploded[[patient]]


    inf.transf = x$fit$transfer[[patient]]
    adj_matrix = DataFrameToMatrix(inf.transf)


    # Augment labels
    lbify = function(w) {
      if (w == 'GL')
        return(w)
      else
        paste(w, ' [', occ.table[w, 'Clonal'], ', ', occ.table[w, 'SubClonal']  , ']', sep = '')
    }

    edgeify = function(w) {
      w[1] = strsplit(w[1], split = ' ')[[1]][1]
      w[2] = strsplit(w[2], split = ' ')[[1]][1]

      edge.table$tab[edge.table$tab$from == w[1] &
                       edge.table$tab$to == w[2], 'Count']
    }

    colnames(adj_matrix) = sapply(colnames(adj_matrix), lbify)
    rownames(adj_matrix) = sapply(rownames(adj_matrix), lbify)

    # print(adj_matrix)

    G = igraph::graph_from_adjacency_matrix(adj_matrix)


    # colors for the nodes...
    igraph::V(G)$color  = "white"

    edLabel = apply(igraph::as_edgelist(G), 1, edgeify)
    edLabel = as.matrix(edLabel, ncol = 1)

    # we work on the layout as well, as before
    lay = igraph::layout.reingold.tilford(G, root = 'GL',  mode = 'all')
    rownames(lay) =  igraph::V(G)$name

    # print(adj_matrix)

    # This topological sort can work only if G has no loops
    if(igraph::is_dag(G))
    {
      TS = wrapTS(adj_matrix)
      for (node in TS)
        lay = fixLayer(adj_matrix, node, lay, offset = 2)
      }

    plot(
      G,
      vertex.frame.color = 'white',
      edge.arrow.size = .5,
      edge.color = 'black',
      edge.label = edLabel,
      layout = lay
    )

    mychar = ifelse(is.na(file), '\u2192', '->')

    title(bquote(bold(.(patient) ~ ':') ~ "information transfer"),
          sub = paste("[n m]: n clonal, m subclonal (driver)          ", mychar, ": k patients"))

  },
  warning = function(w) {
    cat(red('WARNING'), '| ')
    print(w)
    warning(w)
  },
  error = function(w) {
    cat(red('ERROR'), '| ')
    warning(w)
    print(w)

  },
  finally = {

  })

  mylayout.off(lot)
  invisible(NULL)
}



##### Auxiliary function to extract trajectories
fixLayer = function(adj_matrix, w, L, offset = 2)
{
  # all downstream nodes and the highest in the layout
  r.w = reach(MatrixToDataFrame(adj_matrix), w)

  if (length(r.w) == 0)
    return(L)

  l.w = L[w, 2]
  m = max(L[r.w, 2, drop = F])

  if (l.w <= m)
    L[r.w, 2] =  L[r.w, 2] - offset - abs(l.w - m)

  return(L)
}

wrapTS = function(M) {
  tryCatch({
    TS = igraph::topo_sort(igraph::graph_from_adjacency_matrix(M), mode = 'out')$name
    return(TS)
  },
  warning = function(w) {

  },
  error = function(w) {

  },
  finally = {
    return(colnames(M))
  })
}


#
# quartz()
# par(mfrow = c(1, 3))
# revolver_plt_fit_patient(Breast.fit, patient = 'PD9850')
# revolver_plt_trajectories_patient(Breast.fit, patient = 'PD9850')
# revolver_plt_itransfer_patient(Breast.fit, patient = 'PD9850')
