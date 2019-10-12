

get_features_nonM = function(x, patients = x$patients)
{
  Np = length(patients)
  Nd = x$n$drivers

  # =-=-=-=-=-=-=-=-
  # Get all data that we need for driver calls across patients
  # =-=-=-=-=-=-=-=-
  All_drivers = lapply(patients, function(p) {
    samples = Samples(x, p)

    drivers = Drivers(x, p)
    drivers$mean_CCF = apply(drivers[, samples], 1, mean)

    drivers
  })

  All_drivers = Reduce(bind_rows, All_drivers) %>%
    select(variantID, patientID, is.clonal, mean_CCF)

  # =-=-=-=-=-=-=-=-
  # Get all data that we need for trajectories
  # =-=-=-=-=-=-=-=-
  All_trajectories = lapply(patients, function(p) {
    ITransfer(x, p, rank = 1, type = 'drivers', data = 'fits') %>%
      mutate(patientID = p)
  })

  All_trajectories = Reduce(bind_rows, All_trajectories)  %>%
    mutate(
      trajectory = paste0(from, ' --> ', to)
    )

  return(
    list(
      All_drivers = All_drivers,
      All_trajectories = All_trajectories
    )
  )
}



plot_clusters_alternative_ggplot = function(x,
                         cluster_palette = revolver:::distinct_palette_few,
                         driver_palette = revolver:::distinct_palette_many,
                         cutoff_drivers = 10,
                         cutoff_trajectories = 10,
                         arrow.symbol = "-->",
                         ...)
{
  # features = get_features(x, patients)

  drivers = All_drivers %>%
    group_by(variantID) %>%
    filter(n() > cutoff_drivers) %>%
    ungroup()

  trajectories = All_trajectories %>%
    group_by(trajectory) %>%
    filter(n() > cutoff_trajectories) %>%
    ungroup()

  cluster_assignemnts = Cluster(x, patients) %>%
    arrange(cluster)

  # Orderings for this plot: patients and drivers
  patients_ordering = cluster_assignemnts %>% pull(patientID)

  drivers_ordering = Stats_drivers(x, drivers = unique(drivers$variantID)) %>%
    arrange(N_tot) %>%
    pull(variantID)

  # Ordering of patients
  drivers$patientID = factor(drivers$patientID,
                             levels = patients_ordering)

  trajectories$patientID = factor(trajectories$patientID,
           levels = patients_ordering)

  # Ordering of drivers
  drivers$variantID = factor(drivers$variantID,
                             levels = drivers_ordering)


  trajectories$to = factor(trajectories$to,
                           levels = drivers_ordering)

  trajectories = trajectories %>% arrange(desc(to))

  trajectories$trajectory = factor(trajectories$trajectory,
                           levels = unique(trajectories$trajectory))

  pl_mutations = ggplot(drivers,
         aes(
           x = variantID,
           y = patientID,
           fill = mean_CCF,
           color = is.clonal
         )) +
    geom_tile(aes(width = .8, height = .8), size = .5 * cex) +
    scale_color_manual(values = c(`TRUE` = 'red', `FALSE` = NA)) +
    scale_fill_distiller(palette = 'Blues',
                         breaks = seq(0, 1, 0.2),
                         direction = 1) +
    theme_minimal(base_size = 10 * cex) +
    theme(
      legend.position = 'bottom',
      legend.key.size = unit(3, 'mm'),
      axis.text.x = element_text(angle = 90)
    ) +
    guides(
      fill = guide_colorbar("Mean CCF", barwidth = 6)
    )



  pl_trajectories =
    ggplot(trajectories,
                        aes(
                          x = trajectory,
                          y = patientID,
                          fill = to
                        )) +
    geom_tile(aes(width = .8, height = .8), size = 1 * cex) +
    # scale_color_manual(values = c(`TRUE` = 'red', `FALSE` = NA)) +
    scale_fill_manual(values = driver_palette(length(unique(trajectories$to)))) +
    theme_minimal(base_size = 10 * cex) +
    theme(
      legend.position = 'bottom',
      legend.key.size = unit(3, 'mm'),
      axis.text.x = element_text(angle = 90),
      axis.text.y = element_blank()
    ) +
    guides(
      fill = guide_legend("Driver", nrow = 1)
    ) +
      labs(
        y = ''
      )

  # Cluster strip
  strip_cluster_assignemnts = cluster_assignemnts %>%
    reshape2::melt(id = 'patientID')

  strip_cluster_assignemnts$patientID =
    factor(strip_cluster_assignemnts$patientID, levels = patients_ordering)

  pl_clusters = ggplot(strip_cluster_assignemnts,
         aes(y = patientID, x = variable, fill = value)
         ) +
    geom_tile() +
    scale_fill_manual(values = cluster_palette(x$cluster$fits$K)) +
    theme_minimal(base_size = 10 * cex) +
    theme(
      legend.position = 'bottom',
      legend.key.size = unit(3, 'mm'),
      axis.text.x = element_text(angle = 90),
      axis.text.y = element_blank()
    ) +
    guides(fill = guide_legend("Cluster", nrow = 1)) +
    labs(y = '')

  # ggarrange(
  #   pl_mutations,
  #   pl_trajectories,
  #   nrow = 1,
  #   ncol = 2
  # )

  library(cowplot)

  plot_grid(
    pl_mutations,
    pl_trajectories,
    pl_clusters,
    align = "h",
    nrow = 1,
    rel_widths = c(1, 1, .1)
  )

  # =-=-=-=-=-=-=-=-
  # We plot two features tibbles (in matrix format for pheatmap):
  # - Mean CCF per driver per patient;
  # - Clonal events marks.
  # =-=-=-=-=-=-=-=-

  # Mean CCF
  Matrix_mean_CCF = as.matrix(features$Matrix_mean_CCF %>% select(-patientID))
  rownames(Matrix_mean_CCF) = features$Matrix_mean_CCF %>% pull(patientID)

  # Clonal events, labeled as "x"
  Matrix_clonal_drivers = data.frame(features$Matrix_clonal_drivers %>% select(-patientID),
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
    color = c("white", scols(1:9, "Blues")),
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
