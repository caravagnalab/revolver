
#' @title Compute REVOLVER's evolutionary distance h(P_i, P_j).
#'
#' @details
#' This function can be applied only to an object that contains a fit result,
#' which is then of class \code{"rev_cohort_fit"}.
#' Default values represent the most common configuration that one should use.
#'
#' @param x A \code{"rev_cohort_fit"} object.
#' @param use.GL Whether or not the germline node (GL) should be used to compute the distance.
#' Default is \code{TRUE}.
#' @param transitive.closure Whether or not the transitive closure of the ordering relations
#' should be computed, or not. Default is \code{FALSE}.
#'
#' @return A \code{"rev_cohort_fit"}. object with a new field \code{"cluster"} that contains the computation's results.
#' @export
#' @import crayon
#'
#' @examples
#' data(CRC.cohort)
#' fit = revolver_fit(CRC.cohort)
#' fit =revolver_evo_distance(fit)
revolver_evo_distance = function(x, use.GL = TRUE, transitive.closure = FALSE)
{
  asWeightedMatrix = function(df){
    entries.names = unique(unlist(df[ c(1,2)]))

    M = matrix(0, ncol = length(entries.names), nrow = length(entries.names))
    colnames(M) = rownames(M) = entries.names

    for(i in 1:nrow(df))
      M[df[i, 'from'], df[i, 'to']] = M[df[i, 'from'], df[i, 'to']] + 1
    M
  }

  ###########################  ###########################  ###########################
  ###########################  ###########################  ###########################
  ###########################  ###########################  ###########################

  fit.patients = names(x$fit$phylogenies)

  cat(cyan('* Computing REVOLVER\'s evolutionary distance.\n'),
      yellow('\t        Patients    :'), length(fit.patients), '\n',
      yellow('\t        Comparisons :'), length(fit.patients) * (length(fit.patients) - 1) / 2, '\n',
      yellow('\tTransitive Closures :'), ifelse(transitive.closure, green('YES'), red('NO')), '\n',
      yellow('\t          Germlines :'), ifelse(use.GL, green('YES'), red('NO')),'\n\n')

  transfer = lapply(
    fit.patients,
    information.transfer_exploded_nodes,
    x = x,
    transitive.closure = transitive.closure
  )
  names(transfer) = fit.patients
  transfer.matrix = Reduce(rbind, transfer)

  if(!use.GL) {
    cat(cyan('No GL will reduce from'), nrow(transfer.matrix))

    transfer = lapply(transfer, function(w) w[w$from != 'GL', , drop = FALSE])
    transfer.matrix = Reduce(rbind, transfer)

    cat(cyan(' to'), nrow(transfer.matrix), '\n')
  }

  transfer.matrix = asWeightedMatrix(transfer.matrix)

  distances = matrix(0, nrow = length(fit.patients), ncol = length(fit.patients))
  colnames(distances) = rownames(distances) = fit.patients

  i = 0
  idx = 1:length(fit.patients)
  pb = txtProgressBar(min = 0, max = length(fit.patients) * (length(fit.patients) - 1) / 2, style = 3)
  pb.status = getOption('revolver.progressBar', default = TRUE)

  for(p1 in idx){
    for(p2 in p1:length(fit.patients)){

      # update progress bar
      if(pb.status) setTxtProgressBar(pb, i)
      i = i + 1

      if(p1 == p2) next;

      # print(transfer[[p1]])
      # print(transfer[[p2]])

      distances[p1, p2] = evo_distance(transfer.matrix, transfer[[p1]], transfer[[p2]])
    }
  }

  if(!is.null(x$cluster)) cat(red('Overwriting current distance ...\n'))

  x$cluster = NULL
  x$cluster$distances = t(distances)
  x$cluster$transfer = transfer
  x$cluster$dist.obj  = as.dist(x$cluster$distances, upper = TRUE)
  x$cluster$distances.params = c(use.GL, transitive.closure)
  x$cluster$transfer.matrix = transfer.matrix
  names(x$cluster$distances.params) = c('use.GL', 'transitive.closure')

  # cat(cyan('Computing hierarchical clustering. \n.'))
  return(x)
}


#' @title Compute possible clusterings for standard parameters.
#'
#' @details
#' This requires to have computed the distance via function \code{\link{revolver_evo_distance}}.
#' Compute possible clusterings for standard parameters. This helps
#' to determine and select clusters in an output dendogram. This function uses
#' dendogram cutting strategies available through packages such as \code{cluster} and
#' \code{dynamicTreeCut}.
#'
#' @param x A \code{"rev_cohort_fit"} object for which the evolutionary distance has been computed.
#' @param methods Dendogram methods supported by agnes and cluster packages. Default
#' are all of the following: \code{"average"}, \code{"single"}, \code{"complete"}, \code{"ward"}.
#' @param min.group.size Minimum group size for \code{dynamicTreeCut} functions.
#' @param do.plot \code{TRUE} if you want plots to be saved to \code{file}
#' @param file PDF output file.
#'
#' @return None, this function just produces a report to inspect offline.
#' @export
#' @import crayon
#' @import cluster
#' @import dendextend
#' @import dynamicTreeCut
#'
#' @examples
#' data(CRC.cohort)
#' fit = revolver_fit(CRC.cohort)
#' fit = revolver_evo_distance(fit)
#' revolver_infoclustering(fit) # dumped to disk
revolver_infoclustering = function(x,
                                    methods = c( "average", "single", "complete", "ward"),
                                    min.group.size = 2,
                                    do.plot = FALSE,
                                    file = 'Rplot.pdf')
{
  if(is.null(x$cluster)) stop('Did you compute evo?')

  distance = x$cluster$distances
  params = x$cluster$distances.params
  dist.obj = x$cluster$dist.obj

  # HC -- find best method to agglomerate clusters with AGNES
  if(do.plot) {pdf(file, height = 10, width = 10)}

  cat(cyan('[agnes] Hierarchical clustering: agglomerative coefficients (ACs)\n'))
  stats = infoclustering(dist.obj, methods, do.plot = do.plot)
  print(stats$stats)
  max = stats$max
  cat(cyan('Largest AC:'), yellow(max), '\n')

  if(max == 'ward.D2') max = 'ward'

  hc = agnes(dist.obj, method = max)
  dendogram = as.dendrogram(hc)

  ############### Optimal number of clusters with dynamicTreeCut (automatic suggestion)
  # require(dynamicTreeCut)
  cat(cyan('\n[dynamicTreeCut/dendextend] Optimal number of clusters k.'), yellow('\t with min. group size'), ':', min.group.size, '\n')

  cat(
    yellow('\t      cutreeDynamic'), ':',
    split_dendogram(dendogram, hc, distance, 'cutreeDynamic', min.group = min.group.size, do.plot = do.plot)$k
  )

  cat(
    yellow('\n\t  cutreeDynamicTree'), ':',
    split_dendogram(dendogram, hc, distance, 'cutreeDynamicTree', min.group = min.group.size, do.plot = do.plot)$k
  )

  cat(
    yellow('\n\t       cutreeHybrid'), ':',
    split_dendogram(dendogram, hc, distance, 'cutreeHybrid', min.group = min.group.size, do.plot = do.plot)$k
  )

  cat(
    yellow('\n\tstatic (silhouette)'), ':',
    split_dendogram(dendogram, hc, distance, 'static', min.group = min.group.size, do.plot = do.plot)$k,
    '\n'
  )

  if(do.plot){
    dev.off()
    jamPDF(in.files = file, out.file = file, layout = '2x2')
  }

  invisible()
}



#' #' @title Compute hierarchical clustering for a REVOLVER cohort
#' #'
#' #' @details
#' #' Compute hierarchical clustering for a REVOLVER cohort with fit models.
#' #' You might want to see functions \code{\link{revolver_evo_distance}} and
#' #' \code{\link{revolver_infoclustering}} to first compute REVOLVER's
#' #' evolutionary distance, and inspect basic parameter values.
#' #'
#' #'
#' #' @param x A \code{"rev_cohort_fit"} object for which the evolutionary distance has been computed.
#' #' @param hc.method Method for hierarchial clustering, see \code{\link{revolver_infoclustering}}.
#' #' @param split.method Method to cut the dendogram, see \code{\link{revolver_infoclustering}}.
#' #' @param dendogram.type \code{"rectangle"}, \code{"triangle"} or other formats supported by a \code{dendogram} object.
#' #' @param cutoff.features_annotation Minimum number of observations to include a certain feature
#' #' in the plot. The number is absolute; e.g., with \code{n = 3} we annotate only edges that occurr at least
#' #' in 3 patients
#' #' @param min.group.size Minimum group size for \code{dynamicTreeCut} functions.
#' #' @param plot.clusters Plot clusters, default is \code{TRUE} and will print to PDF.
#' #' @param plot.groups Set this to \code{TRUE} if you want to plot also the resulting groups and their models.
#' #' @param file Output file
#' #' @param ... Extra parameters
#' #'
#' #' @return The input \code{x} with a modified field \code{cluster} with results.
#' #' @export
#' #' @import crayon
#' #'
#' #'
#' #' @examples
#' #' data(CRC.cohort)
#' #' fit = revolver_fit(CRC.cohort)
#' #' fit = revolver_evo_distance(fit)
#' #' fit = revolver_cluster(fit) # dumped also to disk
#' revolver_cluster = function(
#'   x,
#'   hc.method = "ward",
#'   split.method = "cutreeHybrid",
#'   dendogram.type = 'rectangle',
#'   cutoff.features_annotation = 3,
#'   min.group.size = 2,
#'   plot.clusters = TRUE,
#'   plot.groups = TRUE,
#'   file = 'REVOLVER-clustering.pdf',
#'   ...
#' )
#' {
#'   if(is.null(x$cluster)) stop('Did you compute evo?')
#'
#'   if(!plot.clusters) cat(cyan('You did not ask for cluster\'s plots, skipping those.\n'))
#'
#'
#'   distance = x$cluster$distances
#'   params = x$cluster$distances.params
#'   dist.obj = x$cluster$dist.obj
#'
#'   # Compute HC and dendogram
#'   cat(cyan('[agnes] Hierarchical clustering\n'))
#'   hc = cluster::agnes(dist.obj, method = hc.method)
#'   dendogram = stats::as.dendrogram(hc)
#'   dendextend::labels_cex(dendogram) = .5 # decrease cex if there are many elements
#'
#'   cat(cyan('\tmethod :'), hc.method, '\n')
#'   cat(cyan('\t    AC :'), hc$ac)
#'
#'   x$cluster$hc = hc
#'   x$cluster$dendogram = dendogram
#'   x$cluster$hc.method = hc.method
#'
#'   ############### Optimal number of clusters with dynamicTreeCut
#'   cat(cyan('\n[dynamicTreeCut/dendextend] Cutting dendogram\n'))
#'
#'   clusters = split_dendogram(
#'     dendogram, hc, distance,
#'     split.method, min.group = min.group.size, do.plot = FALSE)
#'
#'   x$cluster$clusters = clusters$clusters
#'   x$cluster$k = clusters$k
#'   x$cluster$split.method = split.method
#'   x$cluster$labels.colors = clusters$labels.colors
#'
#'   cat(cyan('\tmethod :'), split.method, '-', ifelse(split.method != 'static', 'from dynamicTreeCut', 'via silhouette scoring'), '\n')
#'   cat(cyan('\t   |g| :'), min.group.size, '\n')
#'   cat(cyan('\t     k :'), clusters$k, '\n')
#'
#'   if(plot.clusters)
#'   {
#'     pdf('Dendogram.pdf', width = 10, height = 10)
#'     plot_dendogram(
#'       hc,
#'       dendogram,
#'       clusters$clusters,
#'       plot.type = dendogram.type,
#'       main = paste("REVOLVER clustering of", x$annotation),
#'       sub = paste(
#'         'Agnes clustering (', hc.method, ' method; AC=', round(hc$ac, 3) ,') \nand ',
#'         split.method,' as dendogram cutting method', sep =''),
#'       colors = x$cluster$labels.colors
#'     )
#'
#'     # bannerplot associated to the custering
#'     cluster::bannerplot(hc, main = 'Banner plot', col = c('gainsboro', 'steelblue'))
#'     dev.off()
#'   }
#'
#'   #################### Features plot -- the most important one?
#'   if(plot.clusters)
#'     revolver_featurePlot(x, width = length(x$patients), height = length(x$patients), file = 'Clonal-table.pdf')
#'
#'   # ##############################
#'   # #### Compare against other clusterings via tanglegram
#'   # ##############################
#'   features = revolver.featureMatrix(x)
#'
#'   if(plot.clusters)
#'   {
#'   cat(cyan('[tanglegram] Comparison against other clusterings: '))
#'
#'   # Occurrence of mutations (binarized from avg. CCFs in features$occurrences)
#'   occurrences.binarized = features$occurrences
#'   occurrences.binarized[occurrences.binarized > 0] = 1
#'
#'   occurrences.binarized = as.matrix(occurrences.binarized)
#'   hc2 = cluster::agnes(stats::dist(occurrences.binarized), method = hc.method)
#'   dendogram2 = stats::as.dendrogram(hc2)
#'   dendextend::labels_cex(dendogram2) = 0.5
#'
#'   # Plotting comparative clusterings via tanglegrams
#'   pdf('tanglegram.pdf', width = 10, height = 10)
#'   dend_list <- dendextend::dendlist(dendogram, dendogram2)
#'
#'   dendextend::tanglegram(dendogram, dendogram2,
#'                highlight_distinct_edges = FALSE, # Turn-off dashed lines
#'                common_subtrees_color_lines = TRUE, # Turn-off line colors
#'                common_subtrees_color_branches = TRUE, # Color common branches
#'                cex_main = 1,
#'                main = paste("Binary occurrences"),
#'                sub = paste("entanglement =", round(dendextend::entanglement(dend_list), 4))
#'   )
#'
#'   plot_dendogram(
#'       hc2,
#'       dendogram2,
#'       x$cluster$clusters,
#'       plot.type = dendogram.type,
#'       main = paste("Dendogram of Binary occurrences"),
#'       sub = '',
#'       colors = x$cluster$labels.colors
#'   )
#'
#'   hc3 = cluster::agnes(dist(features$occurrences.clonal.subclonal), method = hc.method)
#'   dendogram3 = stats::as.dendrogram(hc3)
#'   dendextend::labels_cex(dendogram3) = 0.5
#'   dend_list = dendextend::dendlist(dendogram, dendogram3)
#'
#'   dendextend::tanglegram(dendogram, dendogram3,
#'                highlight_distinct_edges = FALSE, # Turn-off dashed lines
#'                common_subtrees_color_lines = TRUE, # Turn-off line colors
#'                common_subtrees_color_branches = TRUE, # Color common branches
#'                cex_main = 1,
#'                main = paste("Clonal/Subclonal occurrences"),
#'                sub = paste("entanglement =", round(dendextend::entanglement(dend_list), 4))
#'   )
#'
#'   plot_dendogram(
#'       hc3,
#'       dendogram3,
#'       x$cluster$clusters,
#'       plot.type = dendogram.type,
#'       main = paste("Dendogram of Clonal/Subclonal"),
#'       sub = '',
#'       colors = x$cluster$labels.colors
#'   )
#'
#'     dev.off()
#'     jamPDF(in.files = 'tanglegram.pdf', out.file = 'tanglegram.pdf', layout = '2x2')
#'     cat(green('DONE\n'))
#'   }
#'
#'
#'   ##############################
#'   # Plotting the evolutionary distance
#'   # # Annotate each sample with the cluster ID
#'   annotations.samples = data.frame(cluster = clusters$clusters)
#'   annotations.samples = annotations.samples[rownames(features$occurrences), , drop = FALSE]
#'
#'   if(plot.clusters)
#'   {
#'
#'     cat(cyan('Plotting the evolutionary distance: '))
#'
#'     pheatmap::pheatmap(
#'       distance,
#'       main = paste("Evolutionary Distance \n ",
#'                    paste('GL:', params['use.GL'], '\n',
#'                          'Transitive closure:', params['transitive.closure'])),
#'       cluster_rows = as.hclust(hc),
#'       cluster_cols = as.hclust(hc),
#'       color = scols(1:100, "YlOrRd"),
#'       border_color = NA,
#'       show_rownames = T,
#'       show_colnames = T,
#'       # treeheight_row = max(hc$height)/3,
#'       # treeheight_col = max(hc$height)/3,
#'       annotation_row = annotations.samples[, 'cluster', drop = FALSE],
#'       annotation_col = annotations.samples[, 'cluster', drop = FALSE],
#'       annotation_colors = list(cluster = x$cluster$labels.colors),
#'       fontsize_row = 3,
#'       fontsize_col = 3,
#'       # treeheight_col = 70,
#'       # treeheight_row = 70,
#'       cellwidth = 4,
#'       cellheight = 4,
#'       # width = 20,
#'       # height = 20,
#'       file = paste('Edistance.pdf')
#'     )
#'     jamPDF(in.files = c("Dendogram.pdf", 'Edistance.pdf'), out.file = 'Dendogram.pdf', layout = '3x1')
#'     cat(green('DONE\n'))
#'   }
#'
#'
#'   if(plot.groups)
#'   {
#'     cat(cyan('Plotting fit models dividided by group.\n'))
#'
#'     groups = split(annotations.samples, f = annotations.samples$cluster)
#'     transfer = x$cluster$transfer
#'
#'     files = sapply(
#'       1:length(groups),
#'       function(w){
#'
#'         fname = paste('Cluster', names(groups)[w], file, sep = '-')
#'
#'         plot(x, patients = rownames(groups[[w]]),
#'              out.file = fname, plot.stat = FALSE, layout = '1x1', palette = 'Dark2')
#'
#'         tr = Reduce(rbind, transfer[rownames(groups[[w]])])
#'         tr.ids = apply(tr, 1, paste, collapse = '~')
#'         tr.ids = table(tr.ids)
#'         tr.ids = as.data.frame(tr.ids)
#'         rownames(tr.ids) = tr.ids$tr.ids
#'
#'         tr = (unique(tr))
#'         rownames(tr) = apply(tr, 1, paste, collapse = '~')
#'
#'         tr = cbind(tr, tr.ids[rownames(tr), ])
#'         rownames(tr) = tr$tr.ids =  NULL
#'         tr$count = tr$Freq
#'         tr$Freq = tr$Freq / nrow(groups[[w]])
#'         tr = tr[order(tr$Freq, decreasing = T), ]
#'         tr = tr[tr$count > 1, ]
#'
#'         w = min(10, nrow(tr))
#'         if(w == 0) return('')
#'
#'         cat(cyan('\n Counts > 1 table\n'))
#'         print(tr[1:w, ])
#'
#'         pdf("data_output.pdf", height = nrow(tr), width = 8.5)
#'
#'         gridExtra::grid.arrange(gridExtra::tableGrob(tr))
#'         dev.off()
#'
#'         jamPDF(
#'           in.files = c("data_output.pdf", fname),
#'           out.file = fname,
#'           layout = '1x1'
#'         )
#'
#'         fname
#'     })
#'
#'     # jamPDF(in.files = files, out.file = 'groups.pdf', layout = '1x1')
#'   }
#'
#'   # Plot consensus model for each group
#'   if(plot.clusters)
#'   {
#'     groups = x$cluster$clusters
#'     for(g in unique(groups))
#'       revolver_plotrj_consensus(
#'         x,
#'         patients = names(groups[groups == g]), min.cutoff = 3, ML = TRUE, file = paste('Con', g, '.pdf', sep = ''),
#'         annotation = paste("Cluster", g),
#'         col.annotation = x$cluster$labels.colors[as.character(g)])
#'
#'     revolver_plotrj_consensus(x, min.cutoff = 3, ML = TRUE, file = 'ConsensusModel.pdf', width = 10, height = 10, annotation = paste("All cohort"))
#'
#'     jamPDF(
#'       in.files = c("ConsensusModel.pdf", paste('Con', unique(groups), '.pdf', sep = '')),
#'       out.file = 'Consensus.pdf',
#'       layout = '1x1'
#'       # page = 'a4'
#'     )
#'
#'
#'     cat(cyan('Assembling final PDF:'), yellow(file))
#'     #
#'     jamPDF(
#'       in.files = c('Clonal-table.pdf', 'Dendogram.pdf', 'tanglegram.pdf'),
#'       out.file = file,
#'       # page = 'a4',
#'       layout = '1x1'
#'     )
#'   }
#'
#'   # cat(cyan('* jamming PDFs'))
#'   # jamPDF(in.files = c('mainClustering2.pdf', 'Edistance.pdf'), out.file = 'mainClustering2.pdf', layout = '2x2')
#'   # jamPDF(c('mainClustering.pdf', 'mainClustering2.pdf'), out.file = 'mainClustering.pdf', layout = '2x1')
#'   # jamPDF(in.files = c('cutreeDynamic.pdf', 'dendextend.pdf'), out.file = 'dendextend.pdf', layout = '1x2')
#'   #
#'   #
#'   # jamPDF(
#'   #   in.files = c('mainClustering.pdf','Clonal-table.pdf',  'dendextend.pdf',  'tanglegram.pdf', 'AC.pdf' ),
#'   #   out.file = file,
#'   #   layout = '1x1'
#'   # )
#'   # cat(green(' DONE\n'))
#'
#'
#'   return(x)
#' }


#' @title Compute hierarchical clustering for a REVOLVER cohort
#'
#' @details
#' Compute hierarchical clustering for a REVOLVER cohort with fit models.
#' You might want to see functions \code{\link{revolver_evo_distance}} and
#' \code{\link{revolver_infoclustering}} to first compute REVOLVER's
#' evolutionary distance, and inspect basic parameter values and the clusters
#' that they allow to identify.
#'
#' @param x A \code{"rev_cohort_fit"} object for which the evolutionary distance has been computed.
#' @param hc.method Method for hierarchial clustering, see \code{\link{revolver_infoclustering}}.
#' @param split.method Method to cut the dendogram, see \code{\link{revolver_infoclustering}}.
#' @param min.group.size Minimum group size for \code{dynamicTreeCut} functions.
#'
#' @return The input \code{x} with a modified field \code{cluster} with results.
#' @export
#' @import crayon
#'
#' @examples
#' data(CRC.cohort)
#' fit = revolver_fit(CRC.cohort)
#' fit = revolver_evo_distance(fit)
#' fit = revolver_cluster(fit) # dumped also to disk
revolver_cluster = function(
  x,
  hc.method = "ward",
  split.method = "cutreeHybrid",
  min.group.size = 2
)
{
  if(is.null(x$cluster)) stop('Did you compute evo?')

  distance = x$cluster$distances
  params = x$cluster$distances.params
  dist.obj = x$cluster$dist.obj

  # Compute HC and dendogram
  cat(cyan('[agnes] Hierarchical clustering\n'))
  hc = cluster::agnes(dist.obj, method = hc.method)
  dendogram = stats::as.dendrogram(hc)
  dendextend::labels_cex(dendogram) = .5 # decrease cex if there are many elements

  cat(cyan('\tmethod :'), hc.method, '\n')
  cat(cyan('\t    AC :'), hc$ac)

  x$cluster$hc = hc
  x$cluster$dendogram = dendogram
  x$cluster$hc.method = hc.method

  ############### Optimal number of clusters with dynamicTreeCut
  cat(cyan('\n[dynamicTreeCut/dendextend] Cutting dendogram\n'))

  clusters = split_dendogram(
    dendogram, hc, distance,
    split.method, min.group = min.group.size, do.plot = FALSE)

  x$cluster$clusters = clusters$clusters
  x$cluster$k = clusters$k
  x$cluster$split.method = split.method
  x$cluster$labels.colors = clusters$labels.colors

  cat(cyan('\tmethod :'), split.method, '-', ifelse(split.method != 'static', 'from dynamicTreeCut', 'via silhouette scoring'), '\n')
  cat(cyan('\t   |g| :'), min.group.size, '\n')
  cat(cyan('\t     k :'), clusters$k, '\n')

  return(x)
}




