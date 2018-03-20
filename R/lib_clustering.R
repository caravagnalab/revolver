infoclustering = function(dist.obj, methods, do.plot = FALSE){

  # HC -- find best method to agglomerate clusters with AGNES
  # require(cluster)
  names(methods) <- methods

  stats = purrr::map_dbl(
    methods,
    function(w) {
      c = agnes(dist.obj, method = w)
      if(do.plot) plot(c, which.plots = 2, cex = .5)
      return(c$ac)
    })

  # Best HC method has highest AC
  max = names(stats)[which.max(stats)]

  return(list(stats = stats, max = max))
}

split_dendogram = function(dendogram, hc, distance, method, min.group,
                           do.plot = FALSE,
                           plot.type = 'rectangle',
                           palette = 'Set1')
{
  # cutreeDynamic
  if(method == 'cutreeDynamic')
  {
    clusters = cutreeDynamic(as.hclust(hc),
                             minClusterSize = min.group,
                             method = 'tree')

    clusters = clusters[order.dendrogram(dendogram)]
    names(clusters) = hc$order.lab
  }

  if(method == 'cutreeDynamicTree')
  {
    clusters = cutreeDynamicTree(
      as.hclust(hc),
      deepSplit = TRUE,
      minModuleSize = min.group)

    clusters = clusters[order.dendrogram(dendogram)]
    names(clusters) = hc$order.lab
  }

  if(method == 'cutreeHybrid')
  {
    w = capture.output({
      clusters = cutreeHybrid(
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
    dend_k = find_k(dendogram, krange = 1:ceiling(length(hcl$labels)/min.group))

    clusters = dend_k$pamobject$clustering
  }

  nclusters = names(clusters)
  clusters = paste('C', clusters, sep = '')
  names(clusters) = nclusters


  if(is.null(clusters)) stop('Unknown split method?')

  k = max(clusters)
  labels.colors = scols(unique(clusters), palette)


  if(do.plot) plot_dendogram(hc, dendogram, clusters, plot.type = plot.type, main = method, colors = labels.colors)

  return(list(k = k, clusters = clusters, labels.colors = labels.colors))
}



### Plotting functions

plot_dendogram = function(hc, dendogram, clusters, palette = 'Set1', plot.type = 'rectangle', main = 'Dendogram', colors = NA, ...) {
  # require(dendextend)


  tryCatch({
    clusters = clusters[hc$order.lab]
    clusters = as.character(clusters)
    clusters[clusters == '0'] = '0: Not assigned'
    names(clusters) = hc$order.lab


    # cat('\nPlotting\n')
    #print(clusters)

    labels = unique(clusters)

    if(all(is.na(colors))) labels.colors = scols(labels, palette)
    else labels.colors = colors

    labels_cex(dendogram) = .5

    colLab <- function(n, groups) {
      if (is.leaf(n)) {
        a <- attributes(n)
        # find group name
        a.group <- clusters[a$label]
        # retrieve the corresponding color
        attr(n, "nodePar") <- c(a$nodePar, list(lab.col = labels.colors[a.group], lab.bg = "grey50", pch = 20))
        attr(n, "frame.plot") <- TRUE
      }
      n
    }
    dendogram <- dendrapply(dendogram, colLab)

    plot(
      dendogram,
      main = main,
      type = plot.type,
      ...)

    legend('topleft', legend = labels, col = labels.colors, pch = 19, bty = 'n')
  },
  warning = function(w) { cat(red('PLOTTING ERROR -- maybe some method returned no clusters?\n', w)); },
  error = function(w) { cat(red('PLOTTING ERROR -- maybe some method returned no clusters?\n', w));  },
  finally = { }
  )
}

plot_tanglegram = function(x, versus = 'binary', hc.method = 'ward', dendogram, hc, file = NA, width = 10, height = 10) {

  if(versus == 'binary') {
    features = revolver.featureMatrix(x)$occurrences
    features[features > 0] = 1
    hc = agnes(dist(as.matrix(features)), method = hc.method)

    dendogram = as.dendrogram(hc)
    labels_cex(dendogram) = 0.5

    main = 'Binary'
  }

  if(versus == 'clonal-subclonal') {
    features = revolver.featureMatrix(x)$occurrences.clonal.subclonal
    hc = agnes(dist(features), method = hc.method)

    dendogram = as.dendrogram(hc)
    labels_cex(dendogram) = 0.5

    main = 'Clonal/subclonal'
  }

  xdendogram = x$cluster$dendogram
  dend_list <- dendlist(xdendogram, dendogram)

  if(!is.na(file)) pdf(file, width = width, height = height)

  tanglegram(xdendogram, dendogram,
             highlight_distinct_edges = FALSE, # Turn-off dashed lines
             common_subtrees_color_lines = TRUE, # Turn-off line colors
             common_subtrees_color_branches = TRUE, # Color common branches
             cex_main = 1,
             main = main,
             sub = paste("entanglement =", round(entanglement(dend_list), 4))
  )

  if(!is.na(file)) dev.off()

  return(list(hc = hc, dendogram = dendogram))
}
