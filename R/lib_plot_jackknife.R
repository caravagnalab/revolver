#' Plot jackknife statistics.
#'
#' @details
#' It plots a heatmap with the result of the jackknife computations.
#'
#' @param x An object obtained as output of \code{\link{revolver_jackknife}}.
#' @param sorting For edges, default sorting is "heterogeneity", which means that we count how many times we find a
#' certain edge (regardless its value estimated via the jackknife). Any different value will render the
#' ordering alphabetical.
#' @param file Output file, because multiple plots are generated this should always be a file.
#' @param ... Parameters forwarded to the \code{pheatmap} call
#'
#' @return Nothing.
#' @export
#'
#' @examples
#' TODO
revolver_plot_jackknife = function(x, sorting = 'heterogeneity',
                                   cutoff.annotate.numbers = 0.6, cutoff.annotate.numEdges = 3,
                                   file = "REVOLVER-jackknife.pdf", jam.PDFs = TRUE, ...) {

  stopifnot(all(!is.null(x$cluster)))
  stopifnot(all(!is.null(x$jackknife)))
  stopifnot(!is.null(file))
  stopifnot(!is.na(file))

  prtHdr('REVOLVER Jackknife plot',
         "Cutoff on the numbers to annotate in each pheatmap:", cutoff.annotate.numbers,
         "Sorting model for the pheatmaps:", sorting,
         "Cutoff on the entries in the edges boxplot:", cutoff.annotate.numbers,
         "File:", file,
         format = '\t')

  # color.Confidence = c('white', 'gainsboro',
  #                              RColorBrewer::brewer.pal(9, 'YlGnBu'))

  # breaks.Confidence = c(0, 0.001, 0.1, 0.20, 0.3, 0.4, 0.5, 0.6, 0.7, 0.80, 1)


  breaks.Confidence = seq(0, 1, 0.05)
  sample.colors = scols(breaks.Confidence, palette = 'YlGnBu')
  color.Confidence = c('gainsboro', sample.colors)




  ######################################################## PLOT 1 -- resamples
  # cat(crayon::cyan('\n\tVisualizing resamples'))
  #
  # pheatmap::pheatmap(x$jackknife$groups,
  #                    cluster_cols = FALSE,
  #                    cluster_rows = FALSE, border = NA,
  #                    color = c('steelblue', 'gainsboro'),
  #                    main = paste('Jackknife resamples'),
  #                    file = 'resamples.pdf',
  #                    show_colnames = FALSE,
  #                    width = 10,
  #                    height = 10
  #                    )

  ######################################################## PLOT 2a -- co-occurrence
  prtTit('Co-clustering matrix')
  print(tibble::as.tibble(x$jackknife$cluster))

  hc = as.hclust(x$cluster$hc)

  labels.colors = x$cluster$labels.colors
  annotations.samples = data.frame(cluster = x$cluster$clusters)

  M.lbl = round(x$jackknife$cluster, 2)
  M.lbl[M.lbl <= cutoff.annotate.numbers] = ''

  pheatmap::pheatmap(x$jackknife$cluster,
                    main = paste('REVOLVER jackknife statistic: co-clustering\n',
                                 'n =', x$jackknife$params$resamples, '(resamples)',
                                 'with p =', x$jackknife$params$leave.out * 100, 'leave out (%)'),
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
                       display_numbers = M.lbl, number_color = 'white', fontsize_number = 4,
                       cluster_rows = as.hclust(hc),
                       cluster_cols = as.hclust(hc),
                       clustering_method = ifelse(hc.method == 'ward', 'ward.D2', hc.method),
                       file = 'Jackknife-cooccurence.pdf'
    )

  ######################################################## PLOT 2b -- boxplot
  J = x$jackknife$cluster
  groups = names(x$cluster$labels.colors)

  Jm = list()

  for(g in groups) {
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

  prtTit('Boxplot (medians)')
  print(medians)

  pdf('Jackknife-boxplot.pdf')
  boxplot(Jm,
          main = paste("REVOLVER jackknife statistic: co-clustering boxplot"),
          xlab = "Co-clustering probability",
          ylab = "Cluster",
          col = colors,
          border = 'black',
          horizontal = TRUE,
          notch = F,
          names = groups
  )
  dev.off()

  ######################################################## PLOT 3 -- edges
  df = x$jackknife$edges

 # Plotting function
  fun = function(M, title) {
    colors = c('gainsboro',
               RColorBrewer::brewer.pal(6, 'YlGnBu'))

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

    if(sorting == 'heterogeneity') {
      M.cols = M.cols[order(M.cols)]
      M.rows = M.rows[order(M.rows)]

      M = M[, names(M.cols)]
      M.lbl = M.lbl[, names(M.cols)]

      M = M[names(M.rows), ]
      M.lbl = M.lbl[names(M.rows), ]

      ann.row = ann.row[names(M.rows), , drop = FALSE]
      ann.col = ann.col[names(M.cols), , drop = FALSE]
    }

    pheatmap::pheatmap(M,
                       main = title,
                       col = color.Confidence,
                       breaks = breaks.Confidence,
                       legend_breaks = breaks.Confidence,
                       annotation_row = ann.row,
                       annotation_col = ann.col,
                       border = NA, cluster_cols = FALSE, cluster_rows = FALSE,
                       fontsize_row = 6, fontsize_col = 8,
                       cellwidth = 10.5, cellheight = 5.5,
                       display_numbers = M.lbl, number_color = 'white', fontsize_number = 4,
                       file = 'Jackknife-edges.pdf',
                       ...
    )
  }

  M = revolver:::DataFrameToMatrix(df)
  M[TRUE] = 0

  for(i in 1:nrow(df))
    M[df[i, 'from'], df[i, 'to']] = round(df[i, 'count'], 2)

  M = M[order(colnames(M)), ]
  M = M[, order(colnames(M))]


  prtTit('Frequency of edges across resamples')
  print(tibble::as.tibble(x$jackknife$cluster))

  # Plot
  fun(M, title = paste('REVOLVER jackknife statistic: edge frequency across resamples\n',
                       'n =', x$jackknife$params$resamples, '(resamples)',
                       'with p =', x$jackknife$params$leave.out * 100, 'leave out (%)')
      )

  ########################################################
  r = x$jackknife$results

  edge.perRunfreq = NULL
  for(p in 1:length(r)) edge.perRunfreq = rbind(edge.perRunfreq, r[[p]]$features$consensus.explosion)

  mean.vals = lapply(
    split(edge.perRunfreq, f = edge.perRunfreq$edge),
    function(w) mean(w$count))

  mean.vals = mean.vals[order(unlist(mean.vals), decreasing = TRUE)]
  selected = mean.vals[mean.vals > cutoff.annotate.numEdges]

  df = data.frame(edge.perRunfreq)
  df = df[df$edge %in% names(selected), ]

  df$edge = factor(df$edge)
  df$count = as.numeric(df$count)

  prtTit('Number of patients per edge across resamples')
  print(tibble::as.tibble(df))

  require(ggplot2)
  n = length(x$patients)
  # nperc = n * (0.1, 0.25, 0.50, 0.75)

  # pdf(file = 'Jackknife-patient2edges.pdf', height = 0.2 * length(selected), width = 3)
  ggplot(df, aes(y = count, x = reorder(edge, count, FUN=median))) +
    # geom_hline(aes(yintercept = 5), colour = 'red', linetype = "longdash") +
    # geom_hline(aes(yintercept = 15), colour = 'red', linetype = "longdash") +
    # geom_hline(aes(yintercept = 25), colour = 'red', linetype = "longdash") +
    geom_boxplot(alpha = .5) +
    scale_fill_manual(values = 'steelblue') +
    coord_flip() +
    labs(x = 'Trajectory', y = 'Count', title = 'Number of patients\nwith edge', subtitle = paste('Total: n = ', n))
  # dev.off()


  ######################################################## Merge PDFs
  if(jam.PDFs)
    xx = jamPDF(
      in.files = c('Jackknife-cooccurence.pdf', 'Jackknife-boxplot.pdf', 'Jackknife-edges.pdf', 'Jackknife-patient2edges.pdf'),
                   out.file = file, layout = '1x1')


}
