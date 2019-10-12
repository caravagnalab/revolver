#' Compute clusters stability via the jackknife.
#' 
#' @description 
#' 
#' For a set of clusters computed via \code{\link{revolver_cluster}}, you can compute
#' their stability via a jackknife. routine. This funcion runs a kind of bootstrap 
#' routine where subset of patients - a desired number - is removed from the cohort and
#' before re-computing the clusters. In this way, the co-clustering probability of each
#' patient is computed, which leads to a mean clustering stability for each one of the
#' original set of clusters and a frequency for the inference of a particular evolutionary
#' trajectory.
#' 
#' A number of functions are available to plot the results from this jackknife analysis.
#' Note that in general if you require a large number of runs (i.e., resamples), this 
#' computation can take some time. This implementation leverages on the \code{easypar} 
#' package to run in parallel all the re-runs, therefore we suggest to run it on a 
#' multi-core machine to appreciate a speed up in the computations.
#'
#' @param cohort A cohort object where fit and clusters have been computed.
#' @param resamples Number of jackknife samples.
#' @param removal A number in \code{[0,1]} for the percentage of samples to leave out in each jackknife iteration.
#' @param cores.ratio Ratio of cores for parallel execution
#' @param options.fit List of parameters for fitting models. See  \code{\link{revolver_fit}}.
#' @param options.clustering List of parameters for clustering with the germline node GL. See \code{\link{revolver_cluster}}.
#'
#' @return A cohort where a new jackknife field contains result from this analysis
#' @export
#' @import easypar
#'
#' @importFrom utils combn
#'
#' @examples
#' \dontrun{
#'  TODO
#' }
#'
revolver_jackknife = function(cohort,
                              resamples = 100,
                              leave.out = 0.1,
                              cores.ratio = 0.9,
                              options.fit = list(initial.solution = NA, transitive.orderings = FALSE, restarts = 0),
                              options.clustering = list(use.GL = TRUE, transitive.closure = FALSE, min.group.size = 3, hc.method = 'ward', cutoff.features_annotation = 3, split.method = 'cutreeHybrid'),
                              dump.partial = FALSE,
                              overwrite = FALSE
)
{
  if(is.null(cohort$cluster)) stop('You want to first compute clusters?')
  stopifnot(leave.out > 0 & leave.out < 1)
  stopifnot(resamples > 1)

  warning("October 2019 - this is to be implemented with the new easypar.")
  return(cohort)
  
  patients = cohort$patients
  npatients = length(patients)

  hasField = !all(is.null(cohort$jackknife))
  toCompute = (overwrite == TRUE && hasField) || !hasField

  if(toCompute)
    pio::pioHdr('REVOLVER Jackknife',
                toPrint = c(
                  `Resamples` = resamples,
                  `Leave out` = leave.out),
                prefix = '\t -')
  else {
    pio::pioHdr('REVOLVER Jackknife')
    cat(red('\tUsing cached computation\n'))
  }


  r = NULL
  if(toCompute)
  {
    pclusters = setup_parallel(cores.ratio)

    ##################### FIT and CLUSTERING
    pio::pioTit("Computing ... (this might take some time)")

    actual.resamples = resamples

    repeat {
      single.run = jackknife_aux_err_function(cohort, actual.resamples, leave.out, options.fit, options.clustering)
      nsingle.run = length(single.run)

      cat(paste('Collected N =', nsingle.run), 'resamples\n')

      if(actual.resamples == 0) break;

      r = append(r, single.run)
      actual.resamples = actual.resamples - nsingle.run
    }

    pio::pioTit(paste0('Jackknife completed, analyzing results: N =', length(r)))
    stop_parallel(pclusters)

    cohort$jackknife$results = r
    cohort$jackknife$params = list(resamples = resamples, leave.out = leave.out)
  }
  else r = cohort$jackknife$results

  ##################### CO-CLUSTER
  # A function to compute co-occurrences of samples in the clustering's results
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

  pio::pioTit('Co-clustering of input patients (stability matrix, non-normalized)')

  co.clustering = matrix(0, nrow = npatients, ncol = npatients)
  rownames(co.clustering) = colnames(co.clustering) = patients

  for(p in 1:length(r)) co.clustering = coocc(r[[p]]$cluster$clusters, co.clustering)

  print(tibble::as.tibble(co.clustering))

  cohort$jackknife$cluster = co.clustering/cohort$jackknife$params$resamples

  # median Jackknife statistics
  pio::pioTit('Median per cluster')

  # print(co.clustering)

  J = cohort$jackknife$cluster
  groups = names(cohort$cluster$labels.colors)

  Jm = list()

  for(g in groups) {
    members = names(cohort$cluster$clusters[cohort$cluster$clusters == g])

    # print(members)
    # print(g)
    # print(J)
    # print(J[members, members])

    data.boxplot = J[members, members]
    data.boxplot = data.boxplot[upper.tri(data.boxplot)]
    Jm = append(Jm, list(data.boxplot))
  }

  # Order by median
  medians = sapply(Jm, median)
  names(medians) = groups

  print(medians)
  cohort$jackknife$cluster.medians = medians

  ##################### EDGE ESTIMATORS
  pio::pioTit('Edge frequency across resamples')

  all.edges = NULL
  for(p in 1:length(r)) all.edges = c(all.edges, r[[p]]$features$consensus.explosion$edge)

  counts = table(all.edges)
  unique.entries = unique(all.edges)

  all.edges = lapply(unique.entries, function(w){
    k = as.numeric(counts[w])
    w = strsplit(w, '~')[[1]]
    data.frame(from = w[1], to = w[2], count = k, stringsAsFactors = FALSE)
  })
  all.edges = Reduce(rbind, all.edges)
  all.edges$count = all.edges$count/cohort$jackknife$params$resamples
  all.edges = all.edges[order(all.edges$count, decreasing = TRUE), ]

  print(tibble::as.tibble(all.edges))
  cohort$jackknife$edges = all.edges



  return(cohort)
}


# The loop with an error-handling function. That's good to iterate
# untill the required number of resamples has been computed. This
# function assumes that all the parallel stuff has been created beforehand
jackknife_aux_err_function = function(cohort, resamples, leave.out, options.fit, options.clustering) {

  r = foreach(num = 1:resamples, .packages = c("crayon", "igraph", 'matrixStats', 'matrixcalc'), .export = ls(globalenv())) %dopar%
  # for(num in 1:resamples)
  {
    # Error handling is explicit here
    tryCatch(
      {
        # Remove patients in "group"
        subsetted.cohort = revolver_deletePatients(
          cohort,
          sample( # randomized subsampling
            cohort$patients,
            round(length(cohort$patients) * leave.out))
        )

        # Get the list of recurrent drivers, and remove the others
        table = clonal.subclonal.table(subsetted.cohort)
        table = table[table$Counts > 1, ]
        subsetted.cohort = revolver_subsetDrivers(subsetted.cohort, rownames(table))

        # Fit the cohort, compute the distance, and the clusters
        fit = revolver_fit(subsetted.cohort, initial.solution = options.fit$initial.solution, transitive.orderings = options.fit$transitive.orderings, restarts = options.fit$restarts)
        fit = revolver_evo_distance(fit, use.GL = options.clustering$use.GL, transitive.closure = options.clustering$transitive.closure)
        fit = revolver_cluster(fit,
                               # plot.clusters = FALSE, plot.groups = FALSE,
                               hc.method = options.clustering$hc.method,
                               min.group.size = options.clustering$min.group.size,
                               # cutoff.features_annotation = options.clustering$cutoff.features_annotation,
                               split.method  = options.clustering$split.method)

        # Results from the features
        features = revolver.featureMatrix(fit)
        result = list(cluster = fit$cluster, features = features)

        # if(dump.partial) save(result, file = paste('Jackknife-', paste(sample(LETTERS, 10, TRUE), collapse = ''), '.RData', sep = ''))
        result
      },
      error = function(e) { # Errors intercepted output a NULL
        # sink("error.jacknife.log")
        # print(e)
        # sink()
        return(NULL)
      },
      finally = { }
    )


  }

  # print(r)

  # Filter entries that are NULL
  mask = sapply(r, function(w) !all(is.null(w)))

  # print(mask)
  # print(result)

  r[mask]
}



