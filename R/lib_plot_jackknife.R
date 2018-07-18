#' Plot jackknife statistics for co-clustering (heatmap).
#'
#' @details
#' It plots a heatmap with the result of the jackknife computations
#' for the co-clustering probability.
#'
#' @param x An object obtained as output of \code{\link{revolver_jackknife}}.
#' @param cutoff.annotate.numbers Numbers above this value will also be written in the plot.
#' @param file Output file, or NA.
#' @param palette RColorBrewer palette
#' @param ... Parameters forwarded to the \code{pheatmap} call
#'
#' @return None
#' @export
#'
#' @examples
#' data(Breast.fit)
#' revolver_plt_jackknife_coclust(Breast.fit)
revolver_plt_jackknife_coclust = function(x,
                                          cutoff.annotate.numbers = 0.6,
                                          palette = 'YlGnBu',
                                          file = NA,
                                          ...)
{
  obj_has_clusters(x)
  obj_has_jackknife(x)

  stopifnot(palette %in% rownames(RColorBrewer::brewer.pal.info))

  pio::pioHdr(
    'REVOLVER Plot: co-clustering via jackknife (heatmap)',
    toPrint = c(
      `Annotate number above value` = cutoff.annotate.numbers,
      `Output file` = file
    ),
    prefix = '\t'
  )

  breaks.Confidence = seq(0, 1, 0.05)
  sample.colors = scols(breaks.Confidence, palette = palette)
  color.Confidence = c('gainsboro', sample.colors)

  ######################################################## PLOT 2a -- co-occurrence
  pio::pioTit('Co-clustering matrix')
  pio::pioDisp(x$jackknife$cluster)

  hc = as.hclust(x$cluster$hc)

  labels.colors = x$cluster$labels.colors
  annotations.samples = data.frame(cluster = x$cluster$clusters)

  M.lbl = round(x$jackknife$cluster, 2)
  M.lbl[M.lbl <= cutoff.annotate.numbers] = ''

  pheatmap::pheatmap(
    x$jackknife$cluster,
    main = paste(
      'REVOLVER jackknife statistic: co-clustering\n',
      'n =',
      x$jackknife$params$resamples,
      '(resamples)',
      'with p =',
      x$jackknife$params$leave.out * 100,
      'leave out (%)'
    ),
    cellwidth = 10,
    cellheight = 10,
    fontsize_row = 8,
    fontsize_col = 8,
    border = NA,
    col = color.Confidence,
    breaks = breaks.Confidence,
    legend_breaks = breaks.Confidence,
    annotation_row = annotations.samples,
    annotation_col = annotations.samples,
    annotation_colors = list(cluster = labels.colors),
    display_numbers = M.lbl,
    number_color = 'white',
    fontsize_number = 4,
    cluster_rows = as.hclust(hc),
    cluster_cols = as.hclust(hc),
    clustering_method = ifelse(hc.method == 'ward', 'ward.D2', hc.method),
    file = file
  )

  invisible(NULL)
}


#' Plot jackknife statistics for co-clustering (boxplot).
#'
#' @details
#' It plots a boxplot with the result of the jackknife computations
#' for the co-clustering probability.
#'
#' @param x An object obtained as output of \code{\link{revolver_jackknife}}.
#' @param file Output file, or NA.
#'
#' @return None
#' @export
#'
#' @examples
#' data(Breast.fit)
#' revolver_plt_jackknife_coclust_bplot(Breast.fit)
revolver_plt_jackknife_coclust_bplot = function(x, file = NA)
{
  obj_has_clusters(x)
  obj_has_jackknife(x)

  pio::pioHdr(
    'REVOLVER Plot: co-clustering via jackknife (boxplot)',
    toPrint = c(`Output file` = file),
    prefix = '\t'
  )

  ######################################################## PLOT 2b -- boxplot
  J = x$jackknife$cluster
  groups = names(x$cluster$labels.colors)

  Jm = list()

  for (g in groups) {
    members = names(x$cluster$clusters[x$cluster$clusters == g])
    data.boxplot = J[members, members]
    data.boxplot = data.boxplot[upper.tri(data.boxplot)]
    Jm = append(Jm, list(data.boxplot))
  }

  # Order by median
  medians = sapply(Jm, median)
  names(medians) = groups

  medians = x$jackknife$cluster.medians
  Jm = Jm[order(medians)]
  colors = x$cluster$labels.colors[order(medians)]
  groups = groups[order(medians)]

  pio::pioTit('Co-clustering boxplot ')
  pio::pioDisp(medians)

  lot = mylayout.on(file, 1, c(5, 5), cex)

  boxplot(
    Jm,
    main = paste("REVOLVER jackknife statistic: co-clustering boxplot"),
    xlab = "Co-clustering probability",
    ylab = "Cluster",
    col = colors,
    border = 'black',
    horizontal = TRUE,
    notch = F,
    names = groups
  )

  mylayout.off(lot)
  invisible(NULL)
}



#' Plot jackknife statistics for edge detection probability.
#'
#' @details
#' It plots a heatmap with the result of the jackknife computations.
#'
#' @param x An object obtained as output of \code{\link{revolver_jackknife}}.
#' @param sorting Default is "heterogeneity", which means that we count how many times we find a
#' certain edge (regardless its value estimated via the jackknife). Any different value will render the
#' ordering alphabetical.
#' @param file Output file, because multiple plots are generated this should always be a file.
#' @param ... Parameters forwarded to the \code{pheatmap} call
#' @param palette RColorBrewer palette
#'
#' @return Nothing.
#' @export
#'
#' @examples
#' data(Breast.fit)
#' revolver_plt_jackknife_edge_prb(Breast.fit)
revolver_plt_jackknife_edge_prb = function(x,
                                           sorting = 'heterogeneity',
                                           cutoff.annotate.numbers = 0.6,
                                           palette = 'YlGnBu',
                                           file = NA,
                                           ...)
{
  obj_has_clusters(x)
  obj_has_jackknife(x)

  stopifnot(palette %in% rownames(RColorBrewer::brewer.pal.info))

  pio::pioHdr(
    'REVOLVER Plot: edge detection probability via jackknife (heatmap)',
    toPrint = c(
      `Columns ordered by` = sorting,
      `Annotate number above value` = cutoff.annotate.numbers,
      `Output file` = file
    ),
    prefix = '\t'
  )

  breaks.Confidence = seq(0, 1, 0.05)
  sample.colors = scols(breaks.Confidence, palette = palette)
  color.Confidence = c('gainsboro', sample.colors)

  # Plotting function
  fun = function(M, title) {
    colors = c('gainsboro',
               RColorBrewer::brewer.pal(6, palette))

    M.lbl = M
    M.lbl[M.lbl <= cutoff.annotate.numbers] = ''

    M.cols = M
    M.cols[M.cols > 0] = 1
    M.cols = colSums(M.cols)

    M.rows = M
    M.rows[M.rows > 0] = 1
    M.rows = rowSums(M.rows)

    ann.row = data.frame(outgoing = M.rows, row.names = rownames(M))
    ann.col = data.frame(incoming = M.cols, row.names = colnames(M))

    if (sorting == 'heterogeneity') {
      M.cols = M.cols[order(M.cols)]
      M.rows = M.rows[order(M.rows)]

      M = M[, names(M.cols)]
      M.lbl = M.lbl[, names(M.cols)]

      M = M[names(M.rows),]
      M.lbl = M.lbl[names(M.rows),]

      ann.row = ann.row[names(M.rows), , drop = FALSE]
      ann.col = ann.col[names(M.cols), , drop = FALSE]
    }

    pheatmap::pheatmap(
      M,
      main = title,
      col = color.Confidence,
      breaks = breaks.Confidence,
      legend_breaks = breaks.Confidence,
      annotation_row = ann.row,
      annotation_col = ann.col,
      border = NA,
      cluster_cols = FALSE,
      cluster_rows = FALSE,
      fontsize_row = 6,
      fontsize_col = 8,
      cellwidth = 10.5,
      cellheight = 5.5,
      display_numbers = M.lbl,
      number_color = 'white',
      fontsize_number = 4,
      file = file,
      ...
    )
  }

  df = x$jackknife$edges
  M = DataFrameToMatrix(df)
  M[TRUE] = 0

  for (i in 1:nrow(df))
    M[df[i, 'from'], df[i, 'to']] = round(df[i, 'count'], 2)

  M = M[order(colnames(M)),]
  M = M[, order(colnames(M))]

  pio::pioTit('Frequency of edges across resamples')
  pio::pioDisp(x$jackknife$cluster)

  # Plot
  fun(
    M,
    title = paste(
      'REVOLVER jackknife statistic: edge frequency across resamples\n',
      'n =',
      x$jackknife$params$resamples,
      '(resamples)',
      'with p =',
      x$jackknife$params$leave.out * 100,
      'leave out (%)'
    )
  )
}

#' Plot jackknife statistics for the number of patients with an edge.
#'
#' @details
#' It plots a boxplot with the result of the jackknife computations.
#'
#' @param x An object obtained as output of \code{\link{revolver_jackknife}}.
#' @param file Output file, because multiple plots are generated this should always be a file.
#' @param cutoff.annotate.numEdges Minimum mean value to filter some low-freq edges.
#' @param cex Cex graphics.
#'
#' @return None
#' @export
#'
#' @examples
#' data(Breast.fit)
#' revolver_plt_jackknife_edge_counts(Breast.fit)
revolver_plt_jackknife_edge_counts = function(x,
                                              cutoff.annotate.numEdges = 3,
                                              cex = 1,
                                              file = NA)
{
  obj_has_clusters(x)
  obj_has_jackknife(x)

  pio::pioHdr(
    'REVOLVER Plot: edge counts via jackknife (boxplot)',
    toPrint = c(
      `Plot edges with mean count above value` = cutoff.annotate.numEdges,
      `Output file` = file
    ),
    prefix = '\t'
  )

  # Computed on the fly
  r = x$jackknife$results

  edge.perRunfreq = NULL
  for (p in 1:length(r))
    edge.perRunfreq = rbind(edge.perRunfreq, r[[p]]$features$consensus.explosion)

  mean.vals = lapply(split(edge.perRunfreq, f = edge.perRunfreq$edge),
                     function(w)
                       mean(w$count))

  mean.vals = mean.vals[order(unlist(mean.vals), decreasing = TRUE)]
  selected = mean.vals[mean.vals > cutoff.annotate.numEdges]

  df = data.frame(edge.perRunfreq)
  df = df[df$edge %in% names(selected), ]

  df$edge = factor(df$edge)
  df$count = as.numeric(df$count)

  pio::pioTit('Number of patients per edge across resamples')
  pio::pioDisp(df)

  lot = mylayout.on(file, 1, c(5, 5), cex)

  # require(ggplot2)
  n = length(x$patients)

  p = ggplot2::ggplot(df, ggplot2::aes(y = count, x = reorder(edge, count, FUN = median))) +
    ggplot2::geom_boxplot(alpha = .5) +
    ggplot2::scale_fill_manual(values = 'steelblue') +
    ggplot2::coord_flip() +
    ggplot2::labs(
      x = 'Trajectory',
      y = 'Count',
      title = 'Number of patients\nwith edge',
      subtitle = paste('Total: n = ', n)
    )
  print(p)

  mylayout.off(lot)
  invisible(NULL)
}
