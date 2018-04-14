
#' Jackknife estimate
#'
#' @param cohort
#' @param cohort.name
#' @param folder.output
#' @param do.plots
#' @param resamples
#' @param removal
#' @param cores.ratio
#' @param options.fit
#' @param options.clustering
#'
#' @return
#' @export
#' @import parallel
#' @import crayon
#' @import doParallel
#' @import foreach
#'
#' @examples
revolver_jackknife = function(cohort,
                              cohort.name = 'REVOLVER-cohort',
                              folder.output = '.',
                              do.plots = TRUE,
                              resamples = 100,
                              leave.out = 0.1,
                              cores.ratio = 0.9,
                              options.fit = list(initial.solution = NA, transitive.orderings = FALSE, restarts = 10),
                              options.clustering = list(use.GL = TRUE, transitive.closure = FALSE, min.group.size = 3, hc.method = 'ward', cutoff.features_annotation = 3, split.method = 'cutreeHybrid')
)
{
  if(is.null(cohort$cluster)) stop('You want to first compute clusters?')

  patients = cohort$patients
  npatients = length(patients)

  co.clustering = matrix(0, nrow = npatients, ncol = npatients)
  rownames(co.clustering) = colnames(co.clustering) = patients

  groups = lapply(1:resamples, function(r) {
    g = patients[sample(1:npatients, npatients * leave.out)]

    df = data.frame(group = patients, stringsAsFactors = FALSE)
    rownames(df) = patients

    df$group[TRUE] = 1
    df[g, 'group'] = 0

    df$group = as.numeric(df$group)
    df
  })
  groups = Reduce(cbind, groups)
  colnames(groups) = paste("Trial", 1:resamples)

  if(do.plots)
    pheatmap::pheatmap(groups, cluster_cols = FALSE,
                       cluster_rows = FALSE, border = NA,
                       color = c('steelblue', 'gainsboro'),
                       main = 'Resamples',
                       file = paste(cohort.name, '-resamples.pdf', sep = ''),
                       cellwidth = 20,
                       cellheight = 10)

  pclusters = setup_parallel(cores.ratio)

  # r=NULL
  # for(num in 1:resamples)
  r = foreach(num = 1:resamples, .packages = c("crayon", "igraph", 'matrixStats', 'matrixcalc'), .export = ls(globalenv())) %dopar%
  {
    # sink(paste(num, '.txt', sep = ''))

    group = patients[groups[, num] == 0]

    # Remove patients in "group"
    subsetted.cohort = revolver_deletePatients(cohort, group)

    # Get the list of recurrent drivers, and remove the others
    table = revolver:::clonal.subclonal.table(subsetted.cohort)
    table = table[table$Counts > 1, ]
    subsetted.cohort = revolver_subsetDrivers(subsetted.cohort, rownames(table))

    # Fit the cohort, compute the distance, and the clusters
    fit = revolver_fit(subsetted.cohort, initial.solution = options.fit$initial.solution, transitive.orderings = options.fit$transitive.orderings, restarts = options.fit$restarts)

    fit = revolver_evo_distance(fit, use.GL = options.clustering$use.GL, transitive.closure = options.clustering$transitive.closure)

    # cat('evo', file = paste(num,'.txt',sep=''))

    fit = revolver_cluster(fit, plot.clusters = FALSE, plot.groups = FALSE,
                           hc.method = options.clustering$hc.method,
                           min.group.size = options.clustering$min.group.size,
                           cutoff.features_annotation = options.clustering$cutoff.features_annotation,
                           split.method  = options.clustering$split.method)
    # r = append(r, list(fit$cluster))
#
#     print(fit)
#     sink()
    fit$cluster
  }

  cat(cyan('Jackknife completed.'))

  stop_parallel(pclusters)

  coocc = function(l, M) {
    cluster.labels = unique(l)

    for(cl in cluster.labels) {
      cl.assignments = names(l[l == cl])

      pairs = combn(cl.assignments, 2, simplify = F)

      for(p in 1:length(pairs)) {
        M[pairs[[p]][1], pairs[[p]][2]] = M[pairs[[p]][1], pairs[[p]][2]] + 1
        M[pairs[[p]][2], pairs[[p]][1]] = M[pairs[[p]][2], pairs[[p]][1]] + 1
      }
    }
    M
  }

  for(p in 1:length(r)) co.clustering = coocc(r[[p]]$clusters , co.clustering)

  cohort$cluster$jackknife = co.clustering/resamples

  if(do.plots)
  {
    col = RColorBrewer::brewer.pal(8, 'YlGnBu')
    hc = as.hclust(cohort$cluster$hc)
    labels.colors = cohort$cluster$labels.colors
    annotations.samples = data.frame(cluster = cohort$cluster$clusters)

    pheatmap::pheatmap(cohort$cluster$jackknife,
                     main = paste('REVOLVER Jackknife for', cohort.name),
                     cellwidth = 10,
                     cellheight = 10,
                     fontsize_row = 8,
                     fontsize_col = 8,
                     border = NA,
                     col = col,
                     annotation_row = annotations.samples,
                     annotation_col = annotations.samples,
                     annotation_colors = list(cluster = labels.colors),
                     cluster_rows = as.hclust(hc),
                     cluster_cols = as.hclust(hc),
                     clustering_method = ifelse(hc.method == 'ward', 'ward.D2', hc.method),
                     file = paste(cohort.name, '-jackknife.pdf', sep = '')
    )
  }

  # Boxplot of median Jackknife statistics
  J = cohort$cluster$jackknife
  groups = names(cohort$cluster$labels.colors)

  Jm = list()

  for(g in groups) {
    members = names(cohort$cluster$clusters[cohort$cluster$clusters == g])
    data.boxplot = J[members, members]
    data.boxplot = data.boxplot[upper.tri(data.boxplot)]
    Jm = append(Jm, list(data.boxplot))
  }

  # Order by median
  medians = sapply(Jm, median)
  names(medians) = groups

  cat(cyan('Median Jackknife statistics per cluster (stability)'))
  print(medians)
  cohort$cluster$jackknife.medians = medians

  if(do.plots)
  {
    Jm = Jm[order(medians)]
    colors = cohort$cluster$labels.colors[order(medians)]
    groups = groups[order(medians)]

    pdf(paste(cohort.name, '-jackknife-boxplot.pdf', sep = ''))
    boxplot(Jm,
            main = paste("Jackknife statistics with", resamples, 'resamples and', leave.out, '% leave-out', sep =''),
            xlab = "Co-clustering probability (line at 0.75)",
            ylab = "Cluster",
            col = colors,
            border = 'black',
            horizontal = TRUE,
            notch = F,
            names = groups
    )
    abline(v = 0.75, col = 'red', lty = 2, lwd = 2)

    dev.off()
  }

  return(cohort)
}


