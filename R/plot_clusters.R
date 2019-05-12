


#' Title
#'
#' @param x 
#' @param cluster_palette 
#' @param driver_palette 
#' @param CCF_palette 
#' @param cutoff_drivers 
#' @param cutoff_trajectories 
#' @param arrow.symbol 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
plot_clusters = function(x,
                         cluster_palette = revolver:::distinct_palette_few,
                         driver_palette = revolver:::distinct_palette_many,
                         CCF_palette = revolver:::gradient_palette,
                         cutoff_drivers = 10,
                         cutoff_trajectories = 10,
                         arrow.symbol = "-->",
                         ...)
{
  features = get_features(x, x$patients)
  
  # =-=-=-=-=-=-=-=-
  # We plot two features tibbles (in matrix format for pheatmap):
  # - Mean CCF per driver per patient;
  # - Clonal events marks.
  # =-=-=-=-=-=-=-=-
  
  # Mean CCF
  Matrix_mean_CCF = as.matrix(features$Matrix_mean_CCF %>% select(-patientID))
  rownames(Matrix_mean_CCF) = features$Matrix_mean_CCF %>% pull(patientID)
  
  # Clonal events, labeled as "x"
  Matrix_clonal_drivers = data.frame(features$Matrix_clonal_drivers %>% select(-patientID), check.names = F,
                                     stringsAsFactors = FALSE)
  rownames(Matrix_clonal_drivers) = features$Matrix_clonal_drivers %>% pull(patientID)
  
  Matrix_clonal_drivers[Matrix_clonal_drivers == 0] = ''
  Matrix_clonal_drivers[Matrix_clonal_drivers == 1] = 'x'
  
  # Subset matriced according to cutoff_drivers
  drivers_counts = colSums(features$Matrix_drivers %>% select(-patientID))
  drivers_counts = drivers_counts[drivers_counts >= cutoff_drivers]
  
  Matrix_mean_CCF = Matrix_mean_CCF[, names(drivers_counts)]
  Matrix_clonal_drivers = Matrix_clonal_drivers[, names(drivers_counts)]
  
  # For layout purposes, we sort both by driver frequency
  drivers_ordering = order(colSums(Matrix_mean_CCF), # More frequent -> higher CCF --> higher in the ordering
                           decreasing = TRUE)
  
  Matrix_mean_CCF = Matrix_mean_CCF[, drivers_ordering]
  Matrix_clonal_drivers = Matrix_clonal_drivers[, drivers_ordering]
  Matrix_clonal_drivers = Matrix_clonal_drivers[rownames(Matrix_mean_CCF), ]
  
  # =-=-=-=-=-=-=-=-
  # Annotations for clusters
  # =-=-=-=-=-=-=-=-
  
  # Number of clusters
  K = x$cluster$fits$K
  
  # Clustering assignment(s)
  Clusters = data.frame(cluster = x$cluster$fits$labels,
                        stringsAsFactors = FALSE)
  
  # We re-order them alphabetically because all other feature matrices
  # are ordered as such, and then we create a factor so that colors
  # are labelelled consistently by cluster label
  Clusters = Clusters[rownames(Matrix_mean_CCF), , drop = FALSE]
  Clusters$cluster = factor(Clusters$cluster, levels = paste0('C', 1:K))
  
  # Cluster colors
  Cluster_colors = pio:::nmfy(paste0('C', 1:K), cluster_palette(K))
  
  # =-=-=-=-=-=-=-=-
  # Annotations for trajectories
  # =-=-=-=-=-=-=-=-
  Matrix_trajectories = data.frame(
    features$Matrix_trajectories %>% select(-patientID),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rownames(Matrix_trajectories) = features$Matrix_trajectories %>% pull(patientID)
  
  # Subset them according to cutoff_trajectories
  counts_trajectories = colSums(Matrix_trajectories)
  counts_trajectories = counts_trajectories[counts_trajectories >= cutoff_trajectories]
  
  Matrix_trajectories = Matrix_trajectories[, names(counts_trajectories)]
  
  if(nrow(Matrix_clonal_drivers) == 0) stop(
    "Cutoff for drivers", cutoff_drivers, "too stringent, decrease it."
  )
  
  # The colour of the trajectories depends on the destination drivers
  to = strsplit(colnames(Matrix_trajectories), split = ' --> ')
  to = sapply(to, function(w)
    w[2])
  
  to_colors = pio:::nmfy(unique(to),
                         driver_palette(length(unique(to))))
  
  # Actual coloring
  Trajectory_colors = lapply(to,
                             function(w)
                               data.frame(
                                 # `0` = ggplot2::alpha('gainsboro', .2),
                                 `0` = 'white',
                                 `1` = to_colors[w],
                                 stringsAsFactors = FALSE,
                                 check.names = FALSE
                               ))
  names(Trajectory_colors) = colnames(Matrix_trajectories)
  
  # We sort trajectories by frequency as well
  # trajectories_ordering = order(counts_trajectories, decreasing = TRUE)
  # Matrix_trajectories = Matrix_trajectories[, trajectories_ordering]
  order_to = order(to, decreasing = TRUE)
  Matrix_trajectories = Matrix_trajectories[, order_to]
  
  # Adjust the arrow symbol as required by the caller
  colnames(Matrix_trajectories) =
    gsub(
      pattern = '-->',
      replacement = arrow.symbol,
      x = colnames(Matrix_trajectories)
    )
  names(Trajectory_colors) = gsub(
    pattern = '-->',
    replacement = arrow.symbol,
    x = names(Trajectory_colors)
  )
  
  
  # =-=-=-=-=-=-=-=-
  # Pheatmap final layouting
  # =-=-=-=-=-=-=-=-
  
  # We extract the clustering etc which are pre-computed, but we also
  # use them to sort the rows of all matrices consistently
  hc = as.hclust(x$cluster$fits$hc)
  
  # Get dendrogeam ordering
  sample_orderings = hc$labels[hc$order]
  
  # =-=-=-=-=-=-=-=-
  # Annotations assembled all together, with their colors
  # =-=-=-=-=-=-=-=-
  
  All_annotations = Reduce(cbind,
                           list(Clusters, Matrix_trajectories))

  All_colors = append(list(`cluster` = Cluster_colors),
                      Trajectory_colors)
  
  pheatmap(
    Matrix_mean_CCF,
    display_numbers = Matrix_clonal_drivers,
    number_color = 'orange',
    color = c("white", CCF_palette(9)),
    breaks = seq(0, 1.1, 0.1),
    cluster_cols = FALSE,
    cluster_rows = hc,
    clustering_method = ifelse(
      x$cluster$parameters$hc.method == 'ward',
      'ward.D2',
      x$cluster$parameters$hc.method
    ),
    treeheight_row = max(hc$height) / 2,
    cellwidth = 10,
    cellheight = 12,
    na_col = 'gainsboro',
    legend = TRUE,
    annotation_legend = FALSE,
    annotation_row = All_annotations,
    annotation_colors = All_colors,
    main = paste0("REVOLVER clusters for ", x$annotation),
    silent = TRUE
  )
}
#
# pheatmap::pheatmap(
#   # features$occurrences,
#   # cluster_rows = as.hclust(hc),
#   # # clustering_distance_rows = dist.obj,
#   # clustering_method = ifelse(hc.method == 'ward', 'ward.D2', hc.method),
#   # color = c("white", scols(1:9, "Blues")),
#   # breaks = seq(0, 1.1, 0.1),
#   # main = "Input data",
#   annotation_row = annotations.samples,
#   # display_numbers = numbers.matrix,
#   # number_color = 'orange',
#   # legend_breaks = unique(annotations.samples$cluster),
#   annotation_legend = FALSE,
#   # annotation_col = annotations.cols,
#   annotation_colors = append(list(cluster = labels.colors), colors.edges.annotation),
#   # cutree_rows = nGroups,
#   treeheight_row = max(hc$height) / 1.5,
#   cellwidth = 10,
#   cellheight = 12,
#   na_col = 'gainsboro',
#   legend = TRUE
# )
#
# library(dendextend)
#
# ggdend = as.ggdend(x$cluster$fits$dendrogra)
