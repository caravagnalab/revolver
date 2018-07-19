
# Features are organized as follows:
# occurrences: average CCF of driver in a sample
# occurrences.clonal.subclonal : each driver is split into clonal or subclonal, and annotated per sample
# clonal.status: 1 if a driver in a sample is clonal
# edges.noexplosion: edges per sample, not exploded
# edges.explosion: edges per sample, exploded
# counts are:
# consensus.noexplosion: counts for edges.noexplosion
# consensus.explosion: counts for edges.explosion
#' @importFrom stats as.hclust
revolver.featureMatrix = function(x, patients = x$patients)
{
  # 0/1 drivers matrices
  occurrences = matrix(0, ncol = length(x$variantIDs.driver), nrow = length(patients))
  colnames(occurrences) = x$variantIDs.driver
  rownames(occurrences) = patients
  clonal.status = occurrences
  clonal.status[clonal.status == 0] = FALSE

  occurrences.clonal.subclonal = matrix(0, ncol = 2 * length(x$variantIDs.driver), nrow = length(patients))
  colnames(occurrences.clonal.subclonal) = c(paste('Clonal', x$variantIDs.driver), paste('Subclonal', x$variantIDs.driver))
  rownames(occurrences.clonal.subclonal) = patients

  # consensus dfs for edges
  consensus.noexplosion = Reduce(rbind, lapply(x$fit$phylogenies[patients], function(w) w$transfer$drivers))
  consensus.explosion = Reduce(rbind, x$fit$transfer[patients])

  e.noexp = unique(DataFrameToEdges(consensus.noexplosion))
  e.exp = unique(DataFrameToEdges(consensus.explosion))

  # 0/1 edge matrices
  edges.noexplosion = matrix(0, ncol = length(e.noexp), nrow = length(patients))
  colnames(edges.noexplosion) = e.noexp
  rownames(edges.noexplosion) = patients

  edges.explosion = matrix(0, ncol = length(e.exp), nrow = length(patients))
  colnames(edges.explosion) = e.exp
  rownames(edges.explosion) = patients

  # annotate every matrix
  for(pat in patients)
  {
    # this patient
    phylo = x$fit$phylogenies[[pat]]
    data = phylo$dataset
    drivers = data[data$is.driver, 'variantID']

    # edges
    edges.noexplosion[pat, DataFrameToEdges(x$fit$phylogenies[[pat]]$transfer$drivers)] = 1
    edges.explosion[pat, DataFrameToEdges(x$fit$transfer[[pat]])] = 1

    # status of each driver
    for(n in drivers)
    {
      clone = as.character(driver.indexOf(phylo, n))
      ccf = phylo$CCF[clone, phylo$samples_ID]

      occurrences[pat, n] = sum(ccf) / phylo$numRegions
      clonal.status[pat, n] = data[
        data$cluster == clone & data$variantID == n & data$is.driver, 'is.clonal']

      if(data[data$cluster == clone & data$variantID == n & data$is.driver, 'is.clonal'])
        occurrences.clonal.subclonal[pat, paste('Clonal', n)] = 1
      else
        occurrences.clonal.subclonal[pat, paste('Subclonal', n)] = 1
    }
  }

  # Get counts
  consensus.noexplosion$edge = DataFrameToEdges(consensus.noexplosion)
  n = table(consensus.noexplosion$edge)
  cnoexp = consensus.noexplosion[!duplicated(consensus.noexplosion$edge), ]
  cnoexp$count = n[cnoexp$edge]
  cnoexp = cnoexp[order(cnoexp$count, decreasing = TRUE) , , drop = FALSE]

  consensus.explosion$edge = DataFrameToEdges(consensus.explosion)
  n = table(consensus.explosion$edge)
  cexp = consensus.explosion[!duplicated(consensus.explosion$edge), ]
  cexp$count = n[cexp$edge]
  cexp = cexp[order(cexp$count, decreasing = TRUE) , , drop = FALSE]

  return(list(
    occurrences = occurrences,
    occurrences.clonal.subclonal = occurrences.clonal.subclonal,
    clonal.status = clonal.status,
    edges.noexplosion = edges.noexplosion,
    edges.explosion = edges.explosion,
    consensus.noexplosion = cnoexp,
    consensus.explosion = cexp
  ))
}


information.transfer_exploded_nodes = function(x, patient, transitive.closure = FALSE)
{
  f = function(model, var, stopVars)
  {
    aux = function(r)
    {
      # cat(r)
      if(r %in% stopVars) return(r) # stop recursion
      if(is.null(r)) return(NULL) # leaf

      c = children(model, r)

      # recursion, reduction etc.
      return(
        Reduce(union,
               sapply(c, aux))
      )
    }

    return(
      Reduce(union, sapply(children(model, var), aux))
    )
  }

  variables = unlist(lapply(x$fit$substitutions[[patient]], function(w) colnames(w)))
  M = x$fit$exploded[[patient]]

  df = Reduce(
    rbind,
    lapply(c(variables, 'GL'),
           function(v)
             expand.grid(
               from = v,
               to = f(M, v, variables),
               stringsAsFactors = FALSE
             )
    ))

  if(transitive.closure){
    df = DataFrameToMatrix(df)
    df = nem::transitive.closure(df, mat = T, loops = FALSE)

    df = MatrixToDataFrame(df)
  }

  return(df)
}

evo_distance = function(transfer, patient.one, patient.two)
{
  trp1 = patient.one
  trp2 = patient.two

  if(nrow(trp1) > 0){
    rownames(trp1) = DataFrameToEdges(trp1)
    trp1 = cbind(trp1, count = apply(trp1, 1, function(w) transfer[w[1],w[2]]))
  }

  if(nrow(trp2) > 0){
    rownames(trp2) = DataFrameToEdges(trp2)
    trp2 = cbind(trp2, count = apply(trp2, 1, function(w) transfer[w[1],w[2]]))
  }

  if(nrow(trp1) == 0 && nrow(trp2) > 0) return(sum(trp2$count))
  if(nrow(trp1) > 0 && nrow(trp2) == 0) return(sum(trp1$count))
  if(nrow(trp1) == 0 && nrow(trp2) == 0) return(0)

  trp1 = split(trp1, f = trp1$to)
  trp2 = split(trp2, f = trp2$to)

  common.variantsID = intersect(names(trp1), names(trp2))

  trp1.specific = setdiff(names(trp1), common.variantsID)
  trp2.specific = setdiff(names(trp2), common.variantsID)

  dist = sum(unlist(lapply(trp1[trp1.specific], function(x) sum(x$count)))) +
    sum(unlist(lapply(trp2[trp2.specific], function(x) sum(x$count))))

  for(c in common.variantsID) dist = dist + abs(sum(trp1[[c]]$count) - sum(trp2[[c]]$count))

  return(dist)
}


infoclustering = function(dist.obj, methods, do.plot = FALSE){

  # HC -- find best method to agglomerate clusters with AGNES
  # require(cluster)
  names(methods) <- methods

  stats = purrr::map_dbl(
    methods,
    function(w) {
      c = cluster::agnes(dist.obj, method = w)
      if(do.plot) plot(c, which.plots = 2, cex = .5)
      return(c$ac)
    })

  # Best HC method has highest AC
  max = names(stats)[which.max(stats)]

  return(list(stats = stats, max = max))
}

#' @importFrom stats order.dendrogram
#' @importFrom utils capture.output
split_dendogram = function(dendogram, hc, distance, method, min.group,
                           do.plot = FALSE,
                           plot.type = 'rectangle',
                           palette = 'Set1')
{
  # cutreeDynamic
  if(method == 'cutreeDynamic')
  {
    clusters = dynamicTreeCut::cutreeDynamic(as.hclust(hc),
                             minClusterSize = min.group,
                             method = 'tree')

    clusters = clusters[order.dendrogram(dendogram)]
    names(clusters) = hc$order.lab
  }

  if(method == 'cutreeDynamicTree')
  {
    clusters = dynamicTreeCut::cutreeDynamicTree(
      as.hclust(hc),
      deepSplit = TRUE,
      minModuleSize = min.group)

    clusters = clusters[order.dendrogram(dendogram)]
    names(clusters) = hc$order.lab
  }

  if(method == 'cutreeHybrid')
  {
    w = capture.output({
      clusters = dynamicTreeCut::cutreeHybrid(
        as.hclust(hc),
        distance,
        minClusterSize = min.group)$labels
      })

    clusters = clusters[order.dendrogram(dendogram)]
    names(clusters) = hc$order.lab
  }

  if(method == 'static')
  {
    # require(dendextend)

    hcl = as.hclust(hc)
    dend_k = dendextend::find_k(dendogram, krange = 1:ceiling(length(hcl$labels)/min.group))

    clusters = dend_k$pamobject$clustering
  }

  nclusters = names(clusters)
  clusters = paste('C', clusters, sep = '')
  names(clusters) = nclusters


  if(is.null(clusters)) stop('Unknown split method?')

  k = length(unique(clusters))
  labels.colors = scols(unique(clusters), palette)


  if(do.plot) plot_dendogram(hc, dendogram, clusters, plot.type = plot.type, main = method, colors = labels.colors)

  return(list(k = k, clusters = clusters, labels.colors = labels.colors))
}



### Plotting functions
