
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
#' 
#' @import crayon
#' @import cluster
#' @import dendextend
#'
#' @examples
#' data(CRC.cohort)
#' fit = revolver_fit(CRC.cohort)
#' fit = revolver_evo_distance(fit)
#' fit = revolver_cluster(fit) # dumped also to disk
revolver_cluster = function(
  x,
  patients = x$patients,
  hc.method = "ward",
  split.method = "cutreeHybrid",
  min.group.size = 2
)
{
  pio::pioHdr(paste0('REVOLVER Clustering - ', x$annotation))
  
  if(!all(has_fits(x, patients)))
    stop("There are no fits for all the requested patients, aborting.")
  
  # The outputs will be stored in this list
  cluster = list()
  
  # =-=-=-=-=-=-=-=-
  # Compute the pairwise distance among each pair of patient
  # which is REVOLVER's main evolutionary distance
  # =-=-=-=-=-=-=-=-
  pio::pioTit("Computing REVOLVER's evolutionary distance from the Information Transfer")
  
  N = length(patients)
  pio::pioStr("\nPatients :", paste0('N = ', N, ' (', N * (N-1) * 0.5, ' comparisons)' ), suffix = '\n\n')

  # Distance matrix is N x N 
  distances = matrix(0, nrow = N, ncol = N)
  colnames(distances) = rownames(distances) = patients
  
  # External progress bar
  i = 0
  pb = txtProgressBar(min = 0, max = N * (N - 1) / 2, style = 3)
  pb.status = getOption('revolver.progressBar', default = TRUE)
  
  E = x$fit$penalty
  
  # O(n^2) comparisons
  for(p1 in 1:N) {
    for(p2 in p1:N) {
      
      # update progress bar
      if(pb.status) setTxtProgressBar(pb, i)
      i = i + 1
      
      if(p1 == p2) next;
      
      # Both patient have a transfer that we take and augment with E
      p1_IT = ITransfer(x, patients[p1], data = 'fits', type = 'drivers') %>%
        left_join(E, by = c('from', 'to'))
      
      p2_IT = ITransfer(x, patients[p2], data = 'fits', type = 'drivers') %>%
        left_join(E, by = c('from', 'to'))
      
      # The distance is the absolute value of non-shared edges
      distances[patients[p1], patients[p2]] = bind_rows(p1_IT, p2_IT) %>%
        distinct() %>%
        pull(count) %>%
        sum()
    }
  }
   
  # Store these objects inside a distances field - requires transpose
  cluster$distances = list(
    matrix = distances, 
    dist_obj = as.dist(t(distances), upper = TRUE)
  )
  
  # =-=-=-=-=-=-=-=-
  # Compute the hierarchical clusters
  # =-=-=-=-=-=-=-=-
  pio::pioTit("Computing Hierarchical Clustering from the distance")
  
  cluster$parameters = 
    list(
      hc.method = hc.method, 
      split.method = split.method, 
      min.group.size = min.group.size
    )
  
  pio::pioStr("\n  Clustering method", hc.method, suffix = '\n')
  pio::pioStr("      Split method", split.method, suffix = '\n')
  pio::pioStr("Minimum group size", min.group.size, suffix = '\n')
  
  # Agnes for HC computation
  hc = cluster::agnes(cluster$distances$dist_obj, method = hc.method)
    
  # Dendrogram
  dendrogram = as.dendrogram(hc)

  # Store these fits
  cluster$fits =
    list(
      hc = hc,
      dendrogram = as.dendrogram(hc)
    )
  
  # And update them splitting the dendrogram
  cluster$fits$K = split_dendrogram(cluster, split.method, min.group.size)$K
  cluster$fits$labels = split_dendrogram(cluster, split.method, min.group.size)$labels
  

  pio::pioStr("\n\nClusters :", paste0('K = ', cluster$fits$K), suffix = '\n')
  pio::pioStr("Cluster size (n)", '', suffix = '\n')  
  print(
    as_tibble(cluster$fits$labels) %>%
      rename(cluster = value) %>%
      group_by(cluster) %>%
      summarise(n = n()) %>%
      arrange(desc(n))
  )
  
  x$cluster = cluster

  return(x)
}



# Function that splits a dendrogram wrapping calls via other packages
# that are specialized in dendrogram-cutting heuristics
split_dendrogram = function(
  cluster,
  method,
  min.group
  )
{
  # Clustering assignments
  labels = NULL
  
  # Cutting methods from the external package
  if (method == 'cutreeDynamic')
  {
    labels = dynamicTreeCut::cutreeDynamic(
      as.hclust(cluster$fits$hc),
      minClusterSize = min.group,
      method = 'tree'
      )
    
    labels = labels[order.dendrogram(cluster$fits$dendrogram)]
    names(labels) = cluster$fits$hc$order.lab
  }
  
  if (method == 'cutreeDynamicTree')
  {
    labels = dynamicTreeCut::cutreeDynamicTree(
      as.hclust(cluster$fits$hc),
      deepSplit = TRUE,
      minModuleSize = min.group)
    
    labels = labels[order.dendrogram(cluster$fits$dendrogram)]
    names(labels) = cluster$fits$hc$order.lab
  }
  
  if (method == 'cutreeHybrid')
  {
    w = capture.output({
      labels = dynamicTreeCut::cutreeHybrid(
        as.hclust(cluster$fits$hc),
        cluster$distances$matrix,
        minClusterSize = min.group)$labels
    })
    
    labels = labels[order.dendrogram(cluster$fits$dendrogram)]
    names(labels) = cluster$fits$hc$order.lab
  }
  
  if (method == 'static')
  {
    K = ceiling(
      length(
        as.hclust(cluster$fits$hc)$labels
        ) / min.group
      )
    
    dend_k = dendextend::find_k(
      cluster$fits$dendrogram, 
      krange = 1:K
      )
    
    labels = dend_k$pamobject$clustering
  }
  
  # Splits are complete, we format some output now
  # checking for errors etc
  if (is.null(labels)) stop('Unknown split method, aborting.')
  
  K = length(unique(labels))
  if(K == 1) message(
    "Clustering split returned only 1 cluster, does this make sense?"
  )
  
  labels = pio:::nmfy(
    names(labels),
    paste('C', labels, sep = '')
  ) 
  
  
  return(
    list(K = K, labels = labels)
    )
}


