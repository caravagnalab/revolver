##########################################
########################################## Exported plotting functions for clustering's results
##########################################

#' @title Plot the dendrogram for REVOLVER's clusters.
#'
#' @details
#' This function can be applied only to an object that contains a fit result
#' and clusters. The function will plot the dendrogram for the results of
#' hierarchical clustering, and will color labels by group assignment.
#'
#' @param x A cohort object with fits and clusters.
#' @param cex Cex for the plot
#' @param dendogram.type Type of dendrogram (rectangle, triangle, circular)
#' @param file Output file, if NA no file is used
#'
#' @return nothing
#' @export
#'
#' @examples
#' TODO
revolver_plt_rdendogram = function(x,
                                   cex = 1,
                                   dendogram.type = 'rectangle',
                                   file = 'REVOLVER-Clusters-Dendogram.pdf')
{
  obj_has_clusters(x)

  args = pio:::nmfy(c('Dendogram type', 'Output file'), c(dendogram.type, file))
  pio::pioHdr('REVOLVER Plot: Dendrogram of REVOLVER clusters', args, prefix = '\t')

  if (!is.na(file))
    pdf(file, width = 10 * cex, height = 10 * cex)

  plot_dendogram(
    hc,
    x$cluster$dendogram,
    x$cluster$cluster,
    plot.type = dendogram.type,
    main = paste("REVOLVER clustering of", x$annotation),
    sub = paste(
      'Agnes: ',
      x$cluster$hc.method,
      ' method with AC =',
      round(hc$ac, 3) ,
      '\nand ',
      'Dendogram cut: ',
      x$cluster$split.method,
      sep = ''
    ),
    colors = x$cluster$labels.colors
  )

  if (!is.na(file))
    dev.off()
  invisible(NULL)
}

#' @title Bannerplot for REVOLVER's clusters.
#'
#' @details
#' This function can be applied only to an object that contains a fit result
#' and clusters. The function will plot the bannerplot for the results of
#' hierarchical clustering.
#'
#' @param x A cohort object with fits and clusters.
#' @param cex Cex for the plot
#' @param file Output file, if NA no file is used
#'
#' @return nothing
#' @export
#'
#' @examples
#' TODO
revolver_plt_rbannerplot = function(x,
                                    cex = 1,
                                    file = 'REVOLVER-Clusters-Bannerplot.pdf')
{
  obj_has_evodistance(x)

  args = pio:::nmfy(c('Output file'), c(file))
  pio::pioHdr('REVOLVER Plot: Bannerplot of REVOLVER clusters', args, prefix = '\t')

  if (!is.na(file))
    pdf(file, width = 10 * cex, height = 10 * cex)

  # bannerplot associated to the custering
  cluster::bannerplot(x$cluster$hc,
                      main = 'Banner plot',
                      col = c('gainsboro', 'steelblue'))

  if (!is.na(file)) dev.off()
  invisible(NULL)
}

#' @title Plot the features of REVOLVER's clusters, in a heatmap.
#'
#' @details
#' This function can be applied only to an object that contains a fit result
#' and clusters. The function will plot the heatmap of all features that define
#' REVOLVER's clusters. These are the data, as well as the repeated evolutionary
#' trajectories detected by the method. This plot is probably the most common one
#' that you want to inspect to analyze your clusters
#'
#' @param x A cohort object with fits and clusters.
#' @param cex Cex for the plot
#' @param cutoff.features_annotation Minimum number of occurrences that a feature needs to have
#' to include it in the heatmap.
#' @param file Output file, if NA no file is used
#'
#' @return nothing
#' @export
#'
#' @examples
#' TODO
revolver_plt_rclusters = function(x,
                                  cex = 1,
                                  cutoff.features_annotation = 2,
                                  file = 'REVOLVER-Clusters-FeaturesHeatmap.pdf')
{
  obj_has_clusters(x)

  args = pio:::nmfy(c('Cutoff to annotate features (min. observations)', 'Output file'), c(cutoff.features_annotation, file))
  pio::pioHdr('REVOLVER Plot: REVOLVER Cluster table with features table', args, prefix = '\t')

  if (!is.na(file))
    pdf(file, width = 10 * cex, height = 10 * cex)

  #################### Features plot -- the most important one?
  revolver_featurePlot(
    x,
    cutoff.features_annotation = cutoff.features_annotation,
    width = length(x$patients) * cex,
    height = length(x$patients) * cex,
    file = file
  )

  if (!is.na(file))
    dev.off()
  invisible(NULL)
}



#' @title Plot the comparison among REVOLVER's clusters and the one obtained
#' by analyzing the occurence of mutations.
#'
#' @details
#' This function can be applied only to an object that contains a fit result
#' and clusters. The function will plot the dendrogram of the alternative
#' clustering method, which will cluster either binary (type = 'binary') or
#' clonal/ subclonal patterns of occurrence (type = 'clonality'). Besides the
#' dendrogram, a tanglegram against REVOLVER's clusters is also computed.
#'
#' @param x A cohort object with fits and clusters.
#' @param cex Cex for the plot
#' @param dendogram.type Type of dendrogram (rectangle, triangle, circular)
#' @param file Output file, if NA no file is used
#' @param jamPDF If TRUE, the produced figures will be jammed. This makes sense
#' only if you set file not NA.
#' @param jam.layout Layout to jam PDFS, e.g. '2x1'
#'
#' @return nothing
#' @export
#'
#' @examples
#' TODO
revolver_plt_compare_dendograms = function(x,
                                           cex = 1,
                                           type = 'binary',
                                           dendogram.type = 'rectangle',
                                           file = paste('REVOLVER-Compare-Dendogram', type, '.pdf', sep = ''),
                                           jamPDF = FALSE,
                                           jam.layout = '2x1')
{
  obj_has_clusters(x)

  args = pio:::nmfy(c('Comparison against', 'Dendrogram type', 'Output file', 'Jam output PDFs', 'Layout for jam'),
                    c(type, dendogram.type, file, jamPDF, jam.layout))
  pio::pioHdr('REVOLVER Dendrogram Plot', args, prefix = '\t')

  if (!is.na(file))
    pdf(file, width = 10 * cex, height = 10 * cex)

  # ##############################
  # #### Compare against other clusterings via tanglegram
  # ##############################
  pio::pioTit(paste('Comparison against alternative clustering via tanglegram'))

  features = revolver.featureMatrix(x)
  hc.method = x$cluster$hc.method

  # REVOLVER results
  hc = x$cluster$hc
  dendogram = x$cluster$dendogram

  # ALternative
  data = hc2 = NULL

  # Binary occurrences
  if (type == 'binary')
  {
    # Occurrence of mutations (binarized from avg. CCFs in features$occurrences)
    data = features$occurrences
    data[data > 0] = 1

    data = as.matrix(data)
    hc2 = cluster::agnes(stats::dist(occurrences.binarized), method = hc.method)
  }

  # Clonal/ subclonal status
  if (type == 'clonality')
  {
    data = features$occurrences.clonal.subclonal
    hc2 = cluster::agnes(stats::dist(features$occurrences.clonal.subclonal),
                         method = hc.method)
  }

  # Dendogram altearnative
  dendogram2 = stats::as.dendrogram(hc2)
  dendextend::labels_cex(dendogram2) = 0.5

  # Plotting comparative clusterings via tanglegrams
  pdf(file, width = 10 * cex, height = 10 * cex)

  dend_list = dendextend::dendlist(dendogram, dendogram2)
  ent = round(dendextend::entanglement(dend_list), 4)

  lbl = ifelse(type == 'binary',
               "Binary Occurrences",
               "Clonal/ Subclonal Occurrences")

  # Tanglegram
  dendextend::tanglegram(
    dendogram,
    dendogram2,
    highlight_distinct_edges = FALSE,
    # Turn-off dashed lines
    common_subtrees_color_lines = TRUE,
    # Turn-off line colors
    common_subtrees_color_branches = TRUE,
    # Color common branches
    cex_main = 1,
    main = lbl,
    sub = paste("entanglement =", ent)
  )

  pio::pioTit(paste('Dendrogram of alternative clustering '))

  # Dendogram of the alternative, coloured with revolver
  plot_dendogram(
    hc2,
    dendogram2,
    x$cluster$clusters,
    plot.type = dendogram.type,
    main = lbl,
    sub = 'Colors: REVOLVER\'s clusters',
    colors = x$cluster$labels.colors
  )

  if (!is.na(file))
    dev.off()
  if (jam.PDF)
    jamPDF(in.files = file,
           out.file = file,
           layout = jam.layout)

  invisible(NULL)
}

#' @title Plot REVOLVER's evolutionary distance, in a heatmap.
#'
#' @details
#' This function can be applied only to an object that contains a fit result
#' and clusters. The function will plot the heatmap of all the distance values
#' of the distance h(.) applied to the samples. The heatmap is annotated with the
#' clustering's dendrogram, as well with coloured labels.
#'
#' @param x A cohort object with fits and clusters.
#' @param cex Cex for the plot
#' @param file Output file, if NA no file is used
#'
#' @return nothing
#' @export
#'
#' @examples
#' TODO
revolver_plt_evodistance = function(x,
                                    cex = 1,
                                    file = 'REVOLVER-Clusters-EvolutionaryDistance.pdf')
{
  obj_has_clusters(x)

  args = pio:::nmfy(c('Output file'), file)
  pio::pioHdr('REVOLVER Evolutionary distance (h) Plot', args, prefix = '\t')

  distance = x$cluster$distances
  params = x$cluster$distances.params
  hc = x$cluster$hc
  features = revolver.featureMatrix(x)

  # Annotate each sample with the cluster ID
  annotations.samples = data.frame(cluster = x$cluster$cluster$clusters)
  annotations.samples = annotations.samples[rownames(features$occurrences), , drop = FALSE]

  pheatmap::pheatmap(
    distance,
    main = paste(
      "Evolutionary Distance \n ",
      paste('GL:', params['use.GL'], '\n',
            'Transitive closure:', params['transitive.closure'])
    ),
    cluster_rows = as.hclust(hc),
    cluster_cols = as.hclust(hc),
    color = scols(1:100, "YlOrRd"),
    border_color = NA,
    show_rownames = T,
    show_colnames = T,
    annotation_row = annotations.samples[, 'cluster', drop = FALSE],
    annotation_col = annotations.samples[, 'cluster', drop = FALSE],
    annotation_colors = list(cluster = x$cluster$labels.colors),
    fontsize_row = 3 * cex,
    fontsize_col = 3 * cex,
    cellwidth = 4 * cex,
    cellheight = 4 * cex,
    file = file
  )

  invisible(NULL)
}


#' @title Plot the consensus model of the trajectories/
#'
#' @details
#' This function can be applied only to an object that contains a fit result
#' and clusters. For each group and for the overall cohort, this function
#' plots the consensus model of the edges, annotated with counts.
#'
#' @param x A cohort object with fits and clusters.
#' @param cex Cex for the plot.
#' @param cutoff.edges_annotation Edges whose frequency in each group is below
#' this threshold will not be plot.
#' @param file Output file, if NA no file is used
#' @param jamPDF If TRUE, the produced figures will be jammed. This makes sense
#' only if you set file not NA.
#' @param jam.layout Layout to jam PDFS, e.g. '2x1'
#'
#' @return nothing
#' @export
#'
#' @examples
#' TODO
revolver_plt_group_trajectories = function(x,
                                           cex = 1,
                                           cutoff.edges_annotation = 3,
                                           file = 'REVOLVER-Clusters-TrajectoriesConsensus.pdf',
                                           jamPDF = FALSE,
                                           jam.layout = '2x1')
{
  obj_has_clusters(x)

  args = pio:::nmfy(c('Cutoff for edges to be annotated (min. occurrences)', 'Output file'), c(cutoff.edges_annotation, file))
  pio::pioHdr('REVOLVER Plot: Trajectories per cluster', args, prefix = '\t')


  groups = x$cluster$clusters

  if (!is.na(file))
    pdf(file, width = 10 * cex, height = 10 * cex)

  for (g in unique(groups))
    revolver_plotrj_consensus(
      x,
      patients = names(groups[groups == g]),
      min.cutoff = cutoff.edges_annotation,
      ML = TRUE,
      file = NA,
      annotation = paste("Cluster", g),
      col.annotation = x$cluster$labels.colors[as.character(g)]
    )

  revolver_plotrj_consensus(
    x,
    min.cutoff = cutoff.edges_annotation,
    ML = TRUE,
    file = NA,
    annotation = paste("All cohort")
  )

  if (!is.na(file))
    dev.off()
  if (jam.PDF)
    jamPDF(in.files = file,
           out.file = file,
           layout = jam.layout)

  invisible(NULL)
}


#' @title Plot model fit grouped by cluster.
#'
#' @details
#' This function can be applied only to an object that contains a fit result
#' and clusters. For each cluster, this function plots the fit of the models
#' in the cluster. All outputs are directed to separate PDF files.
#'
#' @param x A cohort object with fits and clusters.
#' @param cex Cex for the plot.
#' @param file.prefix Each filename will be prefixed by this string,
#' and annotated also with the cluster ID.
#'
#' @return nothing
#' @export
#'
#' @examples
#' TODO
revolver_plt_fit_by_group = function(x,
                                     cex = 1,
                                     file.prefix = 'REVOLVER-Clusters-TrajectoriesConsensus.pdf')
{
  obj_has_clusters(x)

  args = pio:::nmfy(c('Output file'), c(file))
  pio::pioHdr('REVOLVER Plot: Fits divided by cluster', args, prefix = '\t')


  # Annotate each sample with the cluster ID
  annotations.samples = data.frame(cluster = x$cluster$cluster$clusters)
  annotations.samples = annotations.samples[rownames(features$occurrences), , drop = FALSE]

  groups = split(annotations.samples, f = annotations.samples$cluster)
  transfer = x$cluster$transfer

  files = sapply(1:length(groups),
                 function(w) {
                   fname = paste('Cluster', names(groups)[w], file, sep = '-')

                   pio::pioTit(paste("Plotting", fname))

                   plot(
                     x,
                     patients = rownames(groups[[w]]),
                     out.file = fname,
                     plot.stat = FALSE,
                     layout = '1x1',
                     palette = 'Dark2'
                   )

                   tr = Reduce(rbind, transfer[rownames(groups[[w]])])
                   tr.ids = apply(tr, 1, paste, collapse = '~')
                   tr.ids = table(tr.ids)
                   tr.ids = as.data.frame(tr.ids)
                   rownames(tr.ids) = tr.ids$tr.ids

                   tr = (unique(tr))
                   rownames(tr) = apply(tr, 1, paste, collapse = '~')

                   tr = cbind(tr, tr.ids[rownames(tr),])
                   rownames(tr) = tr$tr.ids =  NULL
                   tr$count = tr$Freq
                   tr$Freq = tr$Freq / nrow(groups[[w]])
                   tr = tr[order(tr$Freq, decreasing = T),]
                   tr = tr[tr$count > 1,]

                   w = min(10, nrow(tr))
                   if (w == 0)
                     return('')

                   cat(cyan('\n Counts > 1 table\n'))
                   print(tr[1:w,])

                   pdf("data_output.pdf",
                       height = nrow(tr),
                       width = 8.5)

                   gridExtra::grid.arrange(gridExtra::tableGrob(tr))
                   dev.off()

                   jamPDF(
                     in.files = c("data_output.pdf", fname),
                     out.file = fname,
                     layout = '1x1'
                   )

                   fname
                 })

  invisible(NULL)
}




# distance = x$cluster$distances
# params = x$cluster$distances.params
# dist.obj = x$cluster$dist.obj
# hc.method = x$cluster$hc.method
#
# hc = x$cluster$hc
# dendogram = x$cluster$dendogram
# clusters = x$cluster$cluster
# k = x$cluster$k
# split.method = x$cluster$split.method
# labels.colors =  x$cluster$labels.colors

##############################
# Plotting the evolutionary distance
# # Annotate each sample with the cluster ID
# annotations.samples = data.frame(cluster = clusters$clusters)
# annotations.samples = annotations.samples[rownames(features$occurrences), , drop = FALSE]


##########################################
########################################## Auxiliary plotting functions for clustering
##########################################

plot_dendogram = function(hc,
                          dendogram,
                          clusters,
                          palette = 'Set1',
                          plot.type = 'rectangle',
                          main = 'Dendogram',
                          colors = NA,
                          ...) {
  # require(dendextend)


  tryCatch({
    clusters = clusters[hc$order.lab]
    clusters = as.character(clusters)
    clusters[clusters == '0'] = '0: Not assigned'
    names(clusters) = hc$order.lab


    # cat('\nPlotting\n')
    #print(clusters)

    labels = unique(clusters)

    if (all(is.na(colors)))
      labels.colors = scols(labels, palette)
    else
      labels.colors = colors

    dendextend::labels_cex(dendogram) = .5

    colLab <- function(n, groups) {
      if (is.leaf(n)) {
        a <- attributes(n)
        # find group name
        a.group <- clusters[a$label]
        # retrieve the corresponding color
        attr(n, "nodePar") <-
          c(a$nodePar,
            list(
              lab.col = labels.colors[a.group],
              lab.bg = "grey50",
              pch = 20
            ))
        attr(n, "frame.plot") <- TRUE
      }
      n
    }
    dendogram <- stats::dendrapply(dendogram, colLab)

    plot(dendogram,
         main = main,
         type = plot.type,
         ...)

    legend(
      'topleft',
      legend = labels,
      col = labels.colors,
      pch = 19,
      bty = 'n'
    )
  },
  warning = function(w) {
    cat(red(
      'PLOTTING ERROR -- maybe some method returned no clusters?\n',
      w
    ))

  },
  error = function(w) {
    cat(red(
      'PLOTTING ERROR -- maybe some method returned no clusters?\n',
      w
    ))

  },
  finally = {

  })
}

plot_tanglegram = function(x, versus = 'binary', hc.method = 'ward', dendogram, hc, file = NA, width = 10, height = 10) {

  if(versus == 'binary') {
    features = revolver.featureMatrix(x)$occurrences
    features[features > 0] = 1
    hc = cluster::agnes(dist(as.matrix(features)), method = hc.method)

    dendogram = stats::as.dendrogram(hc)
    dendextend::labels_cex(dendogram) = 0.5

    main = 'Binary'
  }

  if(versus == 'clonal-subclonal') {
    features = revolver.featureMatrix(x)$occurrences.clonal.subclonal
    hc = cluster::agnes(dist(features), method = hc.method)

    dendogram = stats::as.dendrogram(hc)
    dendextend::labels_cex(dendogram) = 0.5

    main = 'Clonal/subclonal'
  }

  xdendogram = x$cluster$dendogram
  dend_list <- dendextend::dendlist(xdendogram, dendogram)

  if(!is.na(file)) pdf(file, width = width, height = height)

  dendextend::tanglegram(xdendogram, dendogram,
                         highlight_distinct_edges = FALSE, # Turn-off dashed lines
                         common_subtrees_color_lines = TRUE, # Turn-off line colors
                         common_subtrees_color_branches = TRUE, # Color common branches
                         cex_main = 1,
                         main = main,
                         sub = paste("entanglement =", round(dendextend::entanglement(dend_list), 4))
  )

  if(!is.na(file)) dev.off()

  return(list(hc = hc, dendogram = dendogram))
}


########### Plot ML estimates for the information transfer
revolver_plotrj_consensus = function(x, annotation = NA, col.annotation = 'white',
                                     patients = x$patients, min.cutoff = 3, ML = TRUE,
                                     file = NA, ...)
{
  if(is.null(x$fit)) stop('Fit a model first, stopping.')

  features = revolver.featureMatrix(x, patients = patients)
  inf.transf = features$consensus.explosion

  if(min.cutoff >= max(inf.transf$count)) min.cutoff = max(inf.transf$count) - 1
  inf.transf = inf.transf[inf.transf$count > min.cutoff, , drop = FALSE]

  if(ML) {
    inf.transf = split(inf.transf, f = inf.transf$to)
    inf.transf = lapply(inf.transf, function(w) {
      w[which(w$count == max(w$count)), , drop = FALSE]
    })
    inf.transf = Reduce(rbind, inf.transf)
  }

  adj_matrix = DataFrameToMatrix(inf.transf)
  occ.table = clonal.subclonal.table(x)

  # Augment labels
  lbify = function(w){
    if(w == 'GL') return(w)
    else paste(w, ' [', occ.table[w, 'Clonal'], ', ', occ.table[w, 'SubClonal']  ,']', sep = '')
  }

  edgeify = function(w){
    w[1] = strsplit(w[1], split = ' ')[[1]][1]
    w[2] = strsplit(w[2], split = ' ')[[1]][1]

    inf.transf[inf.transf$from == w[1] & inf.transf$to == w[2], 'count']
  }

  colnames(adj_matrix) = sapply(colnames(adj_matrix), lbify)
  rownames(adj_matrix) = sapply(rownames(adj_matrix), lbify)

  G = igraph::graph_from_adjacency_matrix(adj_matrix)

  # colors for the nodes...
  igraph::V(G)$color  = "white"

  edLabel = apply(igraph::as_edgelist(G), 1, edgeify)
  edLabel = as.matrix(edLabel, ncol = 1)

  # we work on the layout -- tree if it has GL
  lay = NULL
  if('GL' %in%  igraph::V(G)$name) {
    lay = igraph::layout.reingold.tilford(G, root = 'GL',  mode = 'all')
    rownames(lay) =  igraph::V(G)$name
  }

  # TS = wrapTS(adj_matrix)
  # for(node in TS) lay = fixLayer(node, lay, offset = 2)

  if(!is.na(file)) pdf(file, ...)
  plot(G,
       vertex.frame.color = 'white',
       edge.arrow.size = .5,
       edge.color = 'black',
       edge.label = edLabel,
       layout = lay)
  title(main = 'Information transfer (consensus)')

  legend('topleft',  title = 'Parameters',
         legend = as.expression(
           c(
             bquote(italic(n) == ~ .(length(patients)) ~ '(patients)' ),
             bquote(rho == ~ .(min.cutoff) ~ '(min. freq.)'),
             bquote(ML == ~ .(as.character(ML)) ~ '(MLE)')
           )),
         pch = 19,
         col = c('forestgreen'),
         bg = add.alpha('forestgreen', .3),
         box.col = 'white'
  )

  legend('topright', legend = annotation, bg = add.alpha(col.annotation, .7), pch = 19,  box.col = 'white')

  if(!is.na(file)) dev.off()
}


#################### Features plot -- the most important one for clustering etc.
revolver_featurePlot = function(x, cutoff.features_annotation = 2,
                                file = NA, width = 20, height = 15, device = 'quartz',
                                annotate.alterations = NA, clinical.covariates = NA)
{

  if(is.null(x$fit)) stop('Fit a model first, stopping.')
  if(is.null(x$cluster)) stop('Cluster your cohort first, stopping.')

  # Obvious things that we need for this plot
  clusters = x$cluster
  features = revolver.featureMatrix(x)
  use.GL = x$cluster$distances.params['use.GL']
  hc = as.hclust(x$cluster$hc)
  hc.method = x$cluster$hc$method
  labels.colors = x$cluster$labels.colors

  if(!all(is.na(annotate.alterations))) {
    features$occurrences = features$occurrences[, annotate.alterations, drop = FALSE]
  }
  features$occurrences = features$occurrences[, names(sort(colSums(features$occurrences), decreasing = T))]

  # Annotate each sample with the cluster ID
  annotations.samples = data.frame(cluster = clusters$clusters)
  annotations.samples = annotations.samples[rownames(features$occurrences), , drop = FALSE]

  clusters = split(annotations.samples, f = annotations.samples$cluster)
  transfer = x$cluster$transfer

  # Featured as edges: get all edges in at least cutoff.features_annotation patients; if required remove GL
  which.features = features$consensus.explosion[features$consensus.explosion$count > cutoff.features_annotation, , drop  = F]

  if(!all(is.na(annotate.alterations))) {
    annotate.alterations = c(annotate.alterations, 'GL') # We add germline

    which.features = which.features[
      which.features$from %in% annotate.alterations &
        which.features$to %in% annotate.alterations, , drop = FALSE]
  }

  if(!use.GL) which.features = which.features[which.features$from != 'GL', , drop = FALSE]
  cat(cyan('Features that will be annotated [ use.GL ='), use.GL, cyan(']'), '\n')
  print(which.features)

  edges.annotation = features$edges.explosion[, which.features$edge, drop = FALSE]
  colnames(edges.annotation) = gsub(pattern = '~', ' \u2192 ', colnames(edges.annotation))

  annotations.samples = cbind(annotations.samples, edges.annotation)

  if(!all(is.na(clinical.covariates))) annotations.samples = cbind(annotations.samples, clinical.covariates)

  # rownames(clinical.covariates) %in% rownames(annotations.samples)

  colors.edges.annotation = sapply(
    colnames(edges.annotation),
    function(w) list(c('0' = 'gainsboro', '1' = 'darkgray')))

  numbers.matrix = features$clonal.status[, colnames(features$occurrences), drop = FALSE]
  numbers.matrix[numbers.matrix == 1] = '\u25a0'
  numbers.matrix[numbers.matrix == 0] = ''


  # device
  if(!is.na(file)) dodev(width = width, height = height, device = device)

  setHook("grid.newpage", function() grid::pushViewport(grid::viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")


  pheatmap::pheatmap(features$occurrences,
                     cluster_cols = F,
                     cluster_rows = as.hclust(hc),
                     # clustering_distance_rows = dist.obj,
                     clustering_method = ifelse(hc.method == 'ward', 'ward.D2', hc.method),
                     color = c("white", scols(1:9, "Blues")),
                     breaks = seq(0, 1.1, 0.1),
                     # main = "Input data",
                     annotation_row = annotations.samples,
                     display_numbers = numbers.matrix,
                     number_color = 'orange',
                     # legend_breaks = unique(annotations.samples$cluster),
                     annotation_legend = FALSE,
                     # annotation_col = annotations.cols,
                     annotation_colors = append(list(cluster = labels.colors), colors.edges.annotation),
                     # cutree_rows = nGroups,
                     treeheight_row = max(hc$height)/1.5,
                     cellwidth = 10, cellheight = 12,
                     na_col = 'gainsboro',
                     legend = TRUE
  )


  setHook("grid.newpage", NULL, "replace")

  params = paste('use.GL =', x$cluster$distances.params['use.GL'], '& transitive.closure =', x$cluster$distances.params['transitive.closure'])

  grid::grid.text(
    bquote(bold('Distance ')~italic(h)~' : '~.(params)~bold('  Clustering : ')~.(hc.method)~ ' / '~.(x$cluster$split.method)~' / k ='~.(x$cluster$k)~""),
    y=-0.07, gp=grid::gpar(fontsize=16))

  grid::grid.text(
    bquote(bold('REVOLVER Clusters : ')~.(x$annotation)),
    y=0.97, gp=grid::gpar(fontsize=16))

  grid::grid.text(
    bquote(bold('Left : ')~.('Evolutionary trajectories')),
    x=-0.07, rot=90, gp=grid::gpar(fontsize=16))

  grid::grid.text(
    bquote(bold('Right : ')~.('Input data (mean CCF value, \u25a0 is clonal)')),
    x=0.97, rot=270, gp=grid::gpar(fontsize=16))

  if(!is.na(file)) udodev(file = file)
}

