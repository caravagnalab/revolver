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


#' @title Compute REVOLVER's evolutionary distance h(P_i, P_j).
#'
#' @details
#' Default values represent the most common configuration that one should use.
#'
#' @param x A 'rev_cohort_fit' object.
#' @param use.GL Whether or not the germline node (GL) should be used to compute the distance.
#' Default is TRUE.
#' @param transitive.closure Whether or not the transitive closure of the ordering relations
#' should be computed, or not. Default is FALSE.
#'
#' @return A 'rev_cohort_fit' object with a new field "cluster" that contains the computation's results.
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

  for(p1 in idx){
    for(p2 in p1:length(fit.patients)){

      # update progress bar
      setTxtProgressBar(pb, i)
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
#' Compute possible clusterings for standard parameters. This helps
#' to determine and select clusters in an output dendogram. This function uses
#' dendogram cutting strategies available through packages such as cluster and
#' dynamicTreeCut.
#'
#' @param x A 'rev_cohort_fit' object for which the evolutionary distance has been computed.
#' @param methods Dendogram methods supported by agnes and cluster packages. Default
#' are c( "average", "single", "complete", "ward").
#' @param min.group.size Minimum group size for dynamicTreeCut functions.
#' @param do.plot TRUE if you want plots to be saved to file "file"
#' @param file PDF output file.
#'
#' @return None
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
  require(dynamicTreeCut)
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



#' @title Compute hierarchical clustering for a REVOLVER cohort with fit models.
#'
#' @details Compute hierarchical clustering for a REVOLVER cohort with fit models.
#'
#' @param x A 'rev_cohort_fit' object for which the evolutionary distance has been computed.
#' @param hc.method Method for hierarchial clustering, see "revolver_infoclustering".
#' @param split.method Method to cut the dendogram, see "revolver_infoclustering".
#' @param dendogram.type "rectangle", "triangle" or other formats supported by a dendogram object.
#' @param cutoff.features_annotation Minimum number of observations to include a certain feature
#' in the plot. The number is absolute; e.g., with n = 3 we annotate only edges that occurr at least
#' in 3 patients
#' @param min.group.size Minimum group size for dynamicTreeCut functions.
#' @param plot.groups Set this to TRUE if you want to plot also the resulting groups and their models.
#' @param file Output file
#' @param ... Extra parameters
#'
#' @return A 'rev_cohort_fit' object mod a new field "cluster" that contains the computation's results.
#' @export
#' @import crayon
#' @import cluster
#' @import dendextend
#' @import grid
#' @import dynamicTreeCut
#' @import pheatmap
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
  dendogram.type = 'rectangle',
  cutoff.features_annotation = 3,
  min.group.size = 2,
  plot.groups = TRUE,
  file = 'REVOLVER-clustering.pdf',
  ...
)
{
  if(is.null(x$cluster)) stop('Did you compute evo?')

  distance = x$cluster$distances
  params = x$cluster$distances.params
  dist.obj = x$cluster$dist.obj

  # Compute HC and dendogram
  cat(cyan('[agnes] Hierarchical clustering\n'))
  hc = agnes(dist.obj, method = hc.method)
  dendogram = as.dendrogram(hc)
  labels_cex(dendogram) = .5 # decrease cex if there are many elements

  cat(cyan('\tmethod :'), hc.method, '\n')
  cat(cyan('\t    AC :'), hc$ac)

  x$cluster$hc = hc
  x$cluster$dendogram = dendogram

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

  pdf('Dendogram.pdf', width = 10, height = 10)
  plot_dendogram(
    hc,
    dendogram,
    clusters$clusters,
    plot.type = dendogram.type,
    main = paste("REVOLVER clustering of", x$annotation),
    sub = paste(
      'Agnes clustering (', hc.method, ' method; AC=', round(hc$ac, 3) ,') \nand ',
      split.method,' as dendogram cutting method', sep =''),
    colors = x$cluster$labels.colors
  )

  # bannerplot associated to the custering
  bannerplot(hc, main = 'Banner plot', col = c('gainsboro', 'steelblue'))
  dev.off()

   #################### Features plot -- the most important one?
  revolver_featurePlot(x, width = length(x$patients), height = length(x$patients), file = 'Clonal-table.pdf')

  # features = revolver.featureMatrix(x)
  # use.GL = x$cluster$distances.params['use.GL']
  #
  # features$occurrences = features$occurrences[, names(sort(colSums(features$occurrences), decreasing = T))]
  #
  # # Annotate each sample with the cluster ID
  # annotations.samples = data.frame(cluster = clusters$clusters)
  # annotations.samples = annotations.samples[rownames(features$occurrences), , drop = FALSE]
  #
  # clusters = split(annotations.samples, f = annotations.samples$cluster)
  # transfer = x$cluster$transfer
  #
  # # Featured as edges: get all edges in at least cutoff.features_annotation patients; if required remove GL
  # which.features = features$consensus.explosion[features$consensus.explosion$count > cutoff.features_annotation, , drop  = F]
  # if(!use.GL) which.features = which.features[which.features$from != 'GL', , drop = FALSE]
  # cat(cyan('Features that will be annotated [ use.GL ='), use.GL, cyan(']'), '\n')
  # print(which.features)
  #
  # edges.annotation = features$edges.explosion[, which.features$edge, drop = FALSE]
  # colnames(edges.annotation) = gsub('~', ' --> ', colnames(edges.annotation))
  # annotations.samples = cbind(annotations.samples, edges.annotation)
  #
  # colors.edges.annotation = sapply(
  #   colnames(edges.annotation),
  #   function(w) list(c('0' = 'gainsboro', '1' = 'darkgray')))
  #
  # numbers.matrix = features$clonal.status[, colnames(features$occurrences), drop = FALSE]
  # numbers.matrix[numbers.matrix == 1] = 'X'
  # numbers.matrix[numbers.matrix == 0] = ''
  #
  # pheatmap::pheatmap(features$occurrences,
  #          cluster_cols = F,
  #          cluster_rows = as.hclust(hc),
  #          # clustering_distance_rows = dist.obj,
  #          clustering_method = ifelse(hc.method == 'ward', 'ward.D2', hc.method),
  #          color = c("white", brewer.pal(9, "Blues")),
  #          breaks = seq(0, 1.1, 0.1),
  #          main = "Clusters of variants' patterns",
  #          annotation_row = annotations.samples,
  #          display_numbers = numbers.matrix,
  #          number_color = 'orange',
  #          # legend_breaks = unique(annotations.samples$cluster),
  #          annotation_legend = FALSE,
  #          # annotation_col = annotations.cols,
  #          annotation_colors = append(
  #            list(cluster = x$clusters$labels.colors),
  #            colors.edges.annotation),
  #          # cutree_rows = nGroups,
  #          treeheight_row = max(hc$height)/1.5,
  #          cellwidth = 10, cellheight = 12,
  #          na_col = 'gainsboro',
  #          legend = TRUE,
  #          file = paste('Clonal-table.pdf')
  #          )

  # ##############################
  # #### Compare against other clusterings via tanglegram
  # ##############################
  cat(cyan('[tanglegram] Comparison against other clusterings: '))
  features = revolver.featureMatrix(x)

  # Occurrence of mutations (binarized from avg. CCFs in features$occurrences)
  occurrences.binarized = features$occurrences
  occurrences.binarized[occurrences.binarized > 0] = 1

  occurrences.binarized = as.matrix(occurrences.binarized)
  hc2 = agnes(dist(occurrences.binarized), method = hc.method)
  dendogram2 = as.dendrogram(hc2)
  labels_cex(dendogram2) = 0.5

  pdf('tanglegram.pdf', width = 10, height = 10)
  dend_list <- dendlist(dendogram, dendogram2)

  tanglegram(dendogram, dendogram2,
             highlight_distinct_edges = FALSE, # Turn-off dashed lines
             common_subtrees_color_lines = TRUE, # Turn-off line colors
             common_subtrees_color_branches = TRUE, # Color common branches
             cex_main = 1,
             main = paste("Binary occurrences"),
             sub = paste("entanglement =", round(entanglement(dend_list), 4))
  )

  plot_dendogram(
    hc2,
    dendogram2,
    x$cluster$clusters,
    plot.type = dendogram.type,
    main = paste("Dendogram of Binary occurrences"),
    sub = '',
    colors = x$cluster$labels.colors
  )

  hc3 = agnes(dist(features$occurrences.clonal.subclonal), method = hc.method)
  dendogram3 = as.dendrogram(hc3)
  labels_cex(dendogram3) = 0.5
  #
  dend_list <- dendlist(dendogram, dendogram3)

  tanglegram(dendogram, dendogram3,
             highlight_distinct_edges = FALSE, # Turn-off dashed lines
             common_subtrees_color_lines = TRUE, # Turn-off line colors
             common_subtrees_color_branches = TRUE, # Color common branches
             cex_main = 1,
             main = paste("Clonal/Subclonal occurrences"),
             sub = paste("entanglement =", round(entanglement(dend_list), 4))
  )

  plot_dendogram(
    hc3,
    dendogram3,
    x$cluster$clusters,
    plot.type = dendogram.type,
    main = paste("Dendogram of Clonal/Subclonal"),
    sub = '',
    colors = x$cluster$labels.colors
  )

  dev.off()
  jamPDF(in.files = 'tanglegram.pdf', out.file = 'tanglegram.pdf', layout = '2x2')
  cat(green('DONE\n'))

  ##############################
  # Plotting the evolutionary distance
  cat(cyan('Plotting the evolutionary distance: '))

  # # Annotate each sample with the cluster ID
  annotations.samples = data.frame(cluster = clusters$clusters)
  annotations.samples = annotations.samples[rownames(features$occurrences), , drop = FALSE]

  # pdf(file, ...)
  pheatmap::pheatmap(
    distance,
    main = paste("Evolutionary Distance \n ",
                 paste('GL:', params['use.GL'], '\n',
                       'Transitive closure:', params['transitive.closure'])),
    cluster_rows = as.hclust(hc),
    cluster_cols = as.hclust(hc),
    color = colorRampPalette(brewer.pal(8, "YlOrRd"))(100),
    border_color = NA,
    show_rownames = T,
    show_colnames = T,
    # treeheight_row = max(hc$height)/3,
    # treeheight_col = max(hc$height)/3,
    annotation_row = annotations.samples[, 'cluster', drop = FALSE],
    annotation_col = annotations.samples[, 'cluster', drop = FALSE],
    annotation_colors = list(cluster = x$cluster$labels.colors),
    fontsize_row = 3,
    fontsize_col = 3,
    # treeheight_col = 70,
    # treeheight_row = 70,
    cellwidth = 4,
    cellheight = 4,
    # width = 20,
    # height = 20,
    file = paste('Edistance.pdf')
  )
  jamPDF(in.files = c("Dendogram.pdf", 'Edistance.pdf'), out.file = 'Dendogram.pdf', layout = '3x1')
  cat(green('DONE\n'))


  if(plot.groups)
  {
    cat(cyan('Plotting fit models dividided by group.\n'))

    groups = split(annotations.samples, f = annotations.samples$cluster)
    transfer = x$cluster$transfer

    files = sapply(
      1:length(groups),
      function(w){

        fname = paste('Cluster', names(groups)[w], file, sep = '-')

        plot(x, patients = rownames(groups[[w]]),
             out.file = fname, plot.stat = FALSE, layout = '1x1', palette = 'Dark2')

        tr = Reduce(rbind, transfer[rownames(groups[[w]])])
        tr.ids = apply(tr, 1, paste, collapse = '~')
        tr.ids = table(tr.ids)
        tr.ids = as.data.frame(tr.ids)
        rownames(tr.ids) = tr.ids$tr.ids

        tr = (unique(tr))
        rownames(tr) = apply(tr, 1, paste, collapse = '~')

        tr = cbind(tr, tr.ids[rownames(tr), ])
        rownames(tr) = tr$tr.ids =  NULL
        tr$count = tr$Freq
        tr$Freq = tr$Freq / nrow(groups[[w]])
        tr = tr[order(tr$Freq, decreasing = T), ]
        tr = tr[tr$count > 1, ]

        w = min(10, nrow(tr))
        if(w == 0) return('')

        cat(cyan('\n Counts > 1 table\n'))
        print(tr[1:w, ])

        pdf("data_output.pdf", height = nrow(tr), width = 8.5)

        require(gridExtra)
        grid.arrange(tableGrob(tr))
        dev.off()

        jamPDF(
          in.files = c("data_output.pdf", fname),
          out.file = fname,
          layout = '1x1'
        )

        fname
    })

    # jamPDF(in.files = files, out.file = 'groups.pdf', layout = '1x1')
  }



  # pheatmap(features$transfer,
  #          cluster_cols = F,
  #          clustering_distance_rows = as.dist(x$cluster$distances),
  #          color = c("white", brewer.pal(9, "Blues")),
  #          main = "Clusters of variants' patterns",
  #          # annotation_row = annotations,
  #          # annotation_col = annotations.cols,
  #          # annotation_colors = annotation_colors,
  #          # cutree_rows = nGroups,
  #          treeheight_row = 150,
  #          cellwidth = 10, cellheight = 12, legend = F
  # )

  # pheatmap(features$transfer)
  # dev.off()

  # Plot consensus model for each group
  groups = x$cluster$clusters
  for(g in unique(groups))
    revolver_plotrj_consensus(
      x,
      patients = names(groups[groups == g]), min.cutoff = 3, ML = TRUE, file = paste('Con', g, '.pdf', sep = ''),
      annotation = paste("Cluster", g),
      col.annotation = x$cluster$labels.colors[as.character(g)])

  revolver_plotrj_consensus(x, min.cutoff = 3, ML = TRUE, file = 'ConsensusModel.pdf', width = 10, height = 10, annotation = paste("All cohort"))

  jamPDF(
    in.files = c("ConsensusModel.pdf", paste('Con', unique(groups), '.pdf', sep = '')),
    out.file = 'Consensus.pdf',
    layout = '1x1'
    # page = 'a4'
  )


  cat(cyan('Assembling final PDF:'), yellow(file))
  #
  jamPDF(
    in.files = c('Clonal-table.pdf', 'Dendogram.pdf', 'tanglegram.pdf'),
    out.file = file,
    # page = 'a4',
    layout = '1x1'
  )

  # cat(cyan('* jamming PDFs'))
  # jamPDF(in.files = c('mainClustering2.pdf', 'Edistance.pdf'), out.file = 'mainClustering2.pdf', layout = '2x2')
  # jamPDF(c('mainClustering.pdf', 'mainClustering2.pdf'), out.file = 'mainClustering.pdf', layout = '2x1')
  # jamPDF(in.files = c('cutreeDynamic.pdf', 'dendextend.pdf'), out.file = 'dendextend.pdf', layout = '1x2')
  #
  #
  # jamPDF(
  #   in.files = c('mainClustering.pdf','Clonal-table.pdf',  'dendextend.pdf',  'tanglegram.pdf', 'AC.pdf' ),
  #   out.file = file,
  #   layout = '1x1'
  # )
  # cat(green(' DONE\n'))


  return(x)
}

# Features are organized as follows:
# occurrences: average CCF of driver in a sample
# occurrences.clonal.subclonal : each driver is split into clonal or subclonal, and annotated per sample
# clonal.status: 1 if a driver in a sample is clonal
# edges.noexplosion: edges per sample, not exploded
# edges.explosion: edges per sample, exploded
# counts are:
# consensus.noexplosion: counts for edges.noexplosion
# consensus.explosion: counts for edges.explosion
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

  G = graph_from_adjacency_matrix(adj_matrix)

  # colors for the nodes...
  V(G)$color  = "white"

  edLabel = apply(as_edgelist(G), 1, edgeify)
  edLabel = as.matrix(edLabel, ncol = 1)

  # we work on the layout -- tree if it has GL
  lay = NULL
  if('GL' %in%  V(G)$name) {
    lay = layout.reingold.tilford(G, root = 'GL',  mode = 'all')
    rownames(lay) =  V(G)$name
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

  setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")


  pheatmap::pheatmap(features$occurrences,
                     cluster_cols = F,
                     cluster_rows = as.hclust(hc),
                     # clustering_distance_rows = dist.obj,
                     clustering_method = ifelse(hc.method == 'ward', 'ward.D2', hc.method),
                     color = c("white", brewer.pal(9, "Blues")),
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

  grid.text(
    bquote(bold('Distance ')~italic(h)~' : '~.(params)~bold('  Clustering : ')~.(hc.method)~ ' / '~.(x$cluster$split.method)~' / k ='~.(x$cluster$k)~""),
    y=-0.07, gp=gpar(fontsize=16))

  grid.text(
    bquote(bold('REVOLVER Clusters : ')~.(x$annotation)),
    y=0.97, gp=gpar(fontsize=16))

  grid.text(
    bquote(bold('Left : ')~.('Evolutionary trajectories')),
    x=-0.07, rot=90, gp=gpar(fontsize=16))

  grid.text(
    bquote(bold('Right : ')~.('Input data (mean CCF value, \u25a0 is clonal)')),
    x=0.97, rot=270, gp=gpar(fontsize=16))

  if(!is.na(file)) udodev(file = file)
}
