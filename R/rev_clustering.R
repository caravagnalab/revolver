
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
#' @importFrom stats as.dist
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

  obj_has_fit(x)

  args = pio:::nmfy(c('Use GL (germline) to compute the distance', 'Use transitive closures to compute the distance'),
                    c(use.GL, transitive.closure))
  pio::pioHdr('REVOLVER evolutionary distance', args, prefix = '\t')

  fit.patients = names(x$fit$phylogenies)

  pio::pioTit(
    paste(
      "Computing distance with", length(fit.patients), "patients,",
      length(fit.patients) * (length(fit.patients) - 1) / 2,
      "comparisons"))

  transfer = lapply(
    fit.patients,
    information.transfer_exploded_nodes,
    x = x,
    transitive.closure = transitive.closure
  )
  names(transfer) = fit.patients
  transfer.matrix = Reduce(rbind, transfer)

  if(!use.GL)
  {
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


  pio::pioDisp(x$cluster$distances)

  return(x)
}


#' @title Compute possible clusterings for standard parameters.
#'
#' @details
#' This requires to have computed the distance via function \code{\link{revolver_evo_distance}}.
#' Compute possible clusterings for standard parameters. This helps
#' to determine and select clusters in an output dendrogram. This function uses
#' dendrogram cutting strategies available through packages such as \code{cluster} and
#' \code{dynamicTreeCut}.
#'
#' @param x A \code{"rev_cohort_fit"} object for which the evolutionary distance has been computed.
#' @param methods Dendogram methods supported by agnes and cluster packages. Default
#' are all of the following: \code{"average"}, \code{"single"}, \code{"complete"}, \code{"ward"}.
#' @param min.group.size Minimum group size for \code{dynamicTreeCut} functions.
#' @param cex Cex for the plot
#' @param do.plot \code{TRUE} if you want plots to be saved to \code{file}
#' @param file PDF output file.
#'
#' @return None, this function just produces a report to be inspected offline.
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
                                    cex = 1,
                                    do.plot = FALSE,
                                    file = ifelse(!do.plot, NA, "REVOLVER-infoClustering.pdf"))
{
  obj_has_fit(x)
  obj_has_evodistance(x)

  args = pio:::nmfy(c('Hierarchical CLustering (HC) methods tested', 'Minimum size of each cluster', "Export plots", "Output file"),
                    c(paste(methods, collapse = ', '), min.group.size, do.plot, file))
  pio::pioHdr('REVOLVER infoClustering', args, prefix = '\t')

  #######

  distance = x$cluster$distances
  params = x$cluster$distances.params
  dist.obj = x$cluster$dist.obj

  # HC -- find best method to agglomerate clusters with AGNES
  if(do.plot)
    pdf(file, height = 10 * cex, width = 10 * cex)

  pio::pioTit('Hierarchical clustering: agglomerative coefficients (ACs) via agnes')
  stats = infoclustering(dist.obj, methods, do.plot = do.plot)
  print(stats$stats)
  max = stats$max
  cat(cyan('Largest AC:'), yellow(max), '\n')

  if(max == 'ward.D2') max = 'ward'

  hc = agnes(dist.obj, method = max)
  dendrogram = as.dendrogram(hc)

  ############### Optimal number of clusters with dynamicTreeCut (automatic suggestion)
  pio::pioTit(paste('Dendrogram cut: optimal number "k" of clusters via dendextend'))

  cat(
    yellow('\t      cutreeDynamic'), ':',
    split_dendrogram(dendrogram, hc, distance, 'cutreeDynamic', min.group = min.group.size, do.plot = do.plot)$k
  )

  cat(
    yellow('\n\t  cutreeDynamicTree'), ':',
    split_dendrogram(dendrogram, hc, distance, 'cutreeDynamicTree', min.group = min.group.size, do.plot = do.plot)$k
  )

  cat(
    yellow('\n\t       cutreeHybrid'), ':',
    split_dendrogram(dendrogram, hc, distance, 'cutreeHybrid', min.group = min.group.size, do.plot = do.plot)$k
  )

  cat(
    yellow('\n\tstatic (silhouette)'), ':',
    split_dendrogram(dendrogram, hc, distance, 'static', min.group = min.group.size, do.plot = do.plot)$k,
    '\n'
  )

  if(do.plot) {
    dev.off()
  }

  invisible(NULL)
}


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
#' @param split.method Method to cut the dendrogram, see \code{\link{revolver_infoclustering}}.
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
  obj_has_fit(x)
  obj_has_evodistance(x)

  args = pio:::nmfy(c('Hierarchical CLustering (HC) method', 'Dendrogram Splitting method', 'Minimum size of each cluster'),
                    c(hc.method, split.method, min.group.size))
  pio::pioHdr('REVOLVER Clustering', args, prefix = '\t')

  #######

  distance = x$cluster$distances
  params = x$cluster$distances.params
  dist.obj = x$cluster$dist.obj

  # Compute HC and dendrogram
  pio::pioTit('Computing Hierarchical clustering with agnes')
  hc = cluster::agnes(dist.obj, method = hc.method)
  dendrogram = stats::as.dendrogram(hc)
  dendextend::labels_cex(dendrogram) = .5 # decrease cex if there are many elements

  cat(cyan('\tmethod :'), hc.method, '\n')
  cat(cyan('\t    AC :'), hc$ac)

  x$cluster$hc = hc
  x$cluster$dendrogram = dendrogram
  x$cluster$hc.method = hc.method

  ############### Optimal number of clusters with dynamicTreeCut
  pio::pioTit('Cutting dendrogram with dendextend')

  clusters = split_dendrogram(
    dendrogram, hc, distance,
    split.method, min.group = min.group.size, do.plot = FALSE)

  x$cluster$clusters = clusters$clusters
  x$cluster$k = clusters$k
  x$cluster$split.method = split.method
  x$cluster$labels.colors = clusters$labels.colors

  cat(cyan('\tmethod :'), split.method, '-', ifelse(split.method != 'static', 'from dynamicTreeCut', 'via silhouette scoring'), '\n')
  cat(cyan('\t   |g| :'), min.group.size, '\n')
  cat(cyan('\t     k :'), clusters$k, '\n')

  pio::pioTit('Clustering assignment and count')

  cat('Groups:')
  print(table(x$cluster$clusters))

  cat('Assignments:')
  pio::pioDisp(x$cluster$clusters)

  return(x)
}




