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
