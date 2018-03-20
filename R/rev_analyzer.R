

#' @title Wrapper for a full analysis with REVOLVER
#' @details
#' This is convenient wrapper to run all REVOLVER's analysis within one command.
#' It takes as input a "rev_cohort" object and performs the following steps:
#' 1) it computes phylogenetic trees (either from binary or CCF data);
#' 2) it performs model fit;
#' 3) if performs clustering.
#' This function requires as parameters all the parameters that would have been
#' passed to the repective functions, where you to run the analysis in single steps.
#' This function saves relevant RData objects for its main steps.
#'
#' @param cohort A "rev_cohort" object.
#' @param type Type of data: "CCF" or "binary". If NA it assumes that "cohort" contains
#' already the tree models for your patients.
#' @param cohort.name Default is 'REVOLVER-cohort', use this to identify output files that
#' will have this string as prefix.
#' @param folder.output Output folder, default is current one ".".
#' @param do.plots Set it to TRUE if you want also plots as output.
#' @param options.trees List of parameters for building input trees. See "revolver_compute_phylogenies"
#' or "revolver_compute_CLtrees"
#' @param options.fit List of parameters for fitting models. See "revolver_fit".
#' @param options.clustering.withGL List of parameters for clustering with the germline node GL. See "revolver_cluster".
#' @param options.clustering.withoutGL List of parameters for clustering without the germline node GL. See "revolver_cluster".
#'
#' @return none
#' @export
#'
#' @examples
#' data(CRC.cohort)
#' revolver_analyzer(CRC.cohort, type = 'binary', do.plots = FALSE)
revolver_analyzer = function(cohort,
                             type = 'CCF',
                             cohort.name = 'REVOLVER-cohort',
                             folder.output = '.',
                             do.plots = TRUE,
                             options.trees = list(sspace.cutoff = 10000, n.sampling = 5000, store.max = 200, overwrite = FALSE),
                             options.fit = list(initial.solution = NA, transitive.orderings = FALSE, restarts = 10),
                             options.clustering.withGL = list(transitive.closure = FALSE, min.group.size = 3, hc.method = 'ward', cutoff.features_annotation = 3, split.method = 'cutreeHybrid'),
                             options.clustering.withoutGL = list(transitive.closure = FALSE, min.group.size = 3, hc.method = 'ward', cutoff.features_annotation = 3, split.method = 'cutreeDynamic')
                             )
{
  if(is.na(type)) cat(red('Warning: "type" is neither CCF nor binary, we will not create input trees and proceed assuming they are already available in "cohort".\n'))

  if(folder.output != '.'){
    dir.create(folder.output)
    setwd(folder.output)
  }

  #################### Create input trees
  if(!is.na(type) && type == 'CCF') {
    for(patient in cohort$patients)
      cohort = revolver_compute_phylogenies(cohort, patient, options = options.trees)

    save(cohort, file = paste(cohort.name, '.phylogenyTrees.RData', sep = ''))
  }

  if(!is.na(type) && type == 'binary') {
    for(patient in cohort$patients)
      cohort = revolver_compute_CLtrees(cohort, patient, options = options.trees)

    save(cohort, file = paste(cohort.name, '.binaryTrees.RData', sep = ''))
  }

  #################### Fit
  fit = revolver_fit(cohort, initial.solution = options.fit$initial.solution, transitive.orderings = options.fit$transitive.orderings, restarts = options.fit$restarts)
  save(fit, file = paste(cohort.name, '.fit.RData', sep = ''))

  if(do.plots) {
    revolver_penaltyPlot(fit)
    plot(fit, out.file = paste(cohort.name, '.fit.pdf', sep = ''), plot.stat = TRUE, layout = '1x1',  palette = 'Dark2')
  }
  #################### Summary table of edges counts
  if(do.plots) {
    table.orderings = rev_table.orderings(fit, intelli.cutoff = 3)
    table.orderings$intelli

    pdf('Orderings.pdf', width = 15)
    aa = lapply(names(table.orderings$intelli), function(w){
      grid.arrange(
        textGrob(w),
        tableGrob(table.orderings$intelli[[w]]))
    })
    dev.off()
  }

  #################### Clustering with GL
  dir.create('With GL')
  setwd('With GL')
  fit = revolver_evo_distance(fit, use.GL = TRUE, transitive.closure = options.clustering.withGL$transitive.closure)

  revolver_infoclustering(fit, min.group.size = options.clustering.withGL$min.group.size, do.plot = do.plots, file = 'ClusterReport.pdf')

  fit = revolver_cluster(fit, plot.groups = do.plots,
                         hc.method = options.clustering.withGL$hc.method,
                         min.group.size = options.clustering.withGL$min.group.size,
                         cutoff.features_annotation = options.clustering.withGL$cutoff.features_annotation,
                         split.method  = options.clustering.withGL$split.method)

  save(fit, file = paste(cohort.name, '.clustering.RData', sep = ''))
  setwd('..')

  #################### Clustering without GL
  dir.create('Without GL')
  setwd('Without GL')
  fit = revolver_evo_distance(fit, use.GL = FALSE, transitive.closure = FALSE)
  revolver_infoclustering(fit, min.group.size = 2, do.plot = TRUE, file = 'ClusterReport.pdf')

  fit = revolver_cluster(fit, plot.groups = do.plots,
                         hc.method = options.clustering.withoutGL$hc.method,
                         min.group.size = options.clustering.withoutGL$min.group.size,
                         cutoff.features_annotation = options.clustering.withoutGL$cutoff.features_annotation,
                         split.method  = options.clustering.withoutGL$split.method)

  save(fit, file = paste(cohort.name, '.clustering.RData', sep = ''))
  setwd('..')

  #################### Closing
  if(folder.output != '.') setwd('..')
}













































