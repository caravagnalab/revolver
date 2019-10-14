



#' @title Compute/ add mutation trees to a REVOLVER cohort.
#'
#' @details
#' This is the analogous of \code{\link{revolver_compute_phylogenies}}, but for binary data and
#' hence it computes Chow-Liu trees using also Suppes' conditions. Parameters have exactly
#' the same meaning of the ones described in \code{\link{revolver_compute_phylogenies}}.
#'
#' @param x An object of class \code{"rev_cohort"}
#' @param patient A patient in the cohort, for which mutation trees are created
#' @param precomputed.trees See \code{\link{revolver_compute_phylogenies}}
#' @param precomputed.scores See \code{\link{revolver_compute_phylogenies}}
#' @param options See \code{\link{revolver_compute_phylogenies}}.
#' @param verbose output type.
#'
#' @return a modififed object of class \code{"rev_cohort"} with available
#' mutation trees for \code{patient}.
#' @export
#' @import crayon
#'
#' @examples
#' data(CRC.cohort)
#' new.cohort = revolver_compute_CLtrees(CRC.cohort, patient = 'adenoma_1')
#' print(new.cohort$phylogenies$adenoma_1)
revolver_compute_CLtrees = function(
  x,
  patient,
  precomputed.trees = NULL,
  precomputed.scores = NULL,
  options = list(sspace.cutoff = 10000,
                 n.sampling = 5000,
                 overwrite = FALSE,
                 store.max = 100),
  verbose = FALSE
)
{
  Original.dataset = x$dataset

  # Prepare output
  if(is.null(x$phylogenies)) x$phylogenies = NULL

  pio::pioHdr(paste("REVOLVER Construct mutational trees (Chow-Liu) for", patient),
              c(
                `Use precomputed trees` = paste(all(!is.null(precomputed.trees))),
                `Maximum state space size to switch from exhaustive to Montecarlo search` = options$sspace.cutoff,
                `Number of Montecarlo samples, if not exhaustive` = options$n.sampling,
                `Overwrite the tree if it is already available` = options$overwrite,
                `Maximumum number of trees to store, if multiple are available` = options$store.max
              ),
              prefix = '\t'
  )

  if(!is.null(x$phylogenies))
  {
    if(patient %in% names(x$phylogenies))
    {
      if(!as.logical(options['overwrite']))
      {
        cat(red('\nModels already available and overwrite is FALSE -- skipping patient.\n'))
        return(x)
      }
    }
  }
  cat('\n')


  # Phylo creation
  x$dataset = x$dataset[x$dataset$patientID  == patient, ]

  samples = names(x$CCF.parser(x$dataset[1, 'CCF']))

  CCF = sapply(x$dataset$CCF, x$CCF.parser)

  if(length(samples) == 1) {
    CCF = matrix(CCF, ncol = 1)
    colnames(CCF) = samples
  } else CCF = t(CCF)

  CCF = apply(CCF, 2, as.numeric)

  if(nrow(x$dataset) == 1) {
    CCF = matrix(CCF, nrow = nrow(x$dataset), ncol = length(samples))
    colnames(CCF) = samples
  }

  rownames(CCF) = rownames(x$dataset)

  x$dataset = x$dataset[, c('Misc', 'patientID', 'variantID', 'cluster', 'is.driver', 'is.clonal')]
  x$dataset = cbind(x$dataset, CCF)

  if(verbose) print(head(x$dataset))

  clusters = clusters.table(x$dataset, samples)
  nclusters = nrow(clusters)

  TREES = SCORES = NULL
  if(!any(is.null(precomputed.trees))){
    cat(yellow('\nPrecomputed trees given as input ... using them.\n'))

    TREES = precomputed.trees
    SCORES = precomputed.scores
  }
  else
  {
    pio::pioTit("Groups/ Clusters in the data of this patient")
    print(clusters)

    if(nclusters == 1)
    {
      cat(red('\nThis model has 1 node, we cannot compute its score.'))

      M = matrix(0, ncol = 1, nrow = 1)
      colnames(M) = rownames(M) = rownames(clusters)

      TREES = append(TREES, list(M))
      SCORES = c(SCORES, 1)
    }
    else
    {
      pio::pioTit("Computing Suppes' extended poset")

      ################### Generate Suppes poset
      POSET = poset(x = clusters, regions = samples)

      # we transform it in the input for the tree sampler
      POSET = lapply(POSET, function(w)
      {
        n = names(w)
        w = as.numeric(w)/length(w)
        names(w) = n
        w
      })

      POSET.EDGES = lapply(POSET, names)
      POSET.EDGES = lapply(names(POSET.EDGES), function(w){
        expand.grid(from = POSET.EDGES[[w]], to = w, stringsAsFactors = FALSE)
      })
      POSET.EDGES = Reduce(rbind, POSET.EDGES)

      pio::pioDisp(POSET.EDGES)

      # # Sampling is carried out if there are more than 'sspace.cutoff' trees, in that case we
      # # sample 'n.sampling' possible trees. Otherwise all possible trees are generated.
      # cat(cyan('[compute_rev_phylogenies] Generating possible solutions (it might take time) ... \n'))

      pio::pioTit("Computing solutions from poset (this might take some time)")

      TREES = all.possible.trees(
        G = POSET.EDGES,
        W = POSET,
        sspace.cutoff = options['sspace.cutoff'],
        n.sampling = options['n.sampling']
      )

      pio::pioTit("Scoring and ranking solutions (this might take some time)")

      # ################## Ranking trees. A tree is good according to the following factors:
      # # 1) the MI among the variables x and y, if they are connected by an edge x --> y [TODO: consider if we really need MI]
      # # 2) the Multinomial probability of edge x --> y in the trees determined by the CCF
      # # 3) for every edge  x --> y, the number of times that the CCF of x is greater than the CCF of y
      # # 3) for every node  x --> y1 ... yK, the number of times that the CCF of x is greater than the sum of the CCFs of y1 ... yK
      binary.data = binarize(x$dataset, samples)

      # 1) MI from binarized data -- different options, with a control sample which avoids 0log(0)
      # • a=0:maximum likelihood estimator (see entropy.empirical)
      # • a=1/2:Jeffreys’ prior; Krichevsky-Trovimov (1991) entropy estimator
      # • a=1:Laplace’s prior
      # • a=1/length(y):Schurmann-Grassberger (1996) entropy estimator
      # • a=sqrt(sum(y))/length(y):minimax prior
      cat('Computing Mutual Information from data\n')
      MI.table = computeMI.table(binary.data, MI.Bayesian.prior = 0, add.control = TRUE)

      pio::pioDisp(MI.table)

      cat('Computing rank\n')
      RANKED = rankTrees(TREES, MI.table, structural.score = NULL)
      TREES = RANKED$TREES
      SCORES = RANKED$SCORES

      TREES = TREES[SCORES > 0]
      SCORES = SCORES[SCORES > 0]
    }
  }

  if(length(TREES) == 0){
    cat(red('No trees found for this patient -- check data? -- returning original cohort.\n'))

    x$dataset = Original.dataset
    return(x)
  }

  pio::pioTit("Creating rev_phylo objects for REVOLVER (this might take some time)")

  x$phylogenies[[patient]] = create_trees_in_revolver_format(options, TREES, SCORES, patient, x$dataset, samples)

  # Restore data
  x$dataset = Original.dataset

  return(x)
}







#' Subset the drivers in the cohort to match the input list.
#'
#' @param cohort An object of class \code{"rev_cohort"}
#' @param list A vector of driver IDs to subset the data to.
#'
#' @return A REVOLVER cohort with only drivers in \code{list}.
#' @export
#'
#' @examples
#' data(CRC.cohort)
#' revolver_subsetDrivers(CRC.cohort, 'APC')
revolver_subsetDrivers = function(cohort, list)
{
  pio::pioHdr('REVOLVER subsetting drivers according to list',
              toPrint = c(
                `Drivers to keep` = paste(list, collapse = ', ')
              ),
              prefix = '\t -')

  current.drivers = rownames(clonal.subclonal.table(cohort))
  toDelete = setdiff(current.drivers, list)

  cat(cyan('Removing:'), paste(toDelete, collapse = ', '))

  cohort$dataset[cohort$dataset$variantID %in% toDelete, 'is.driver'] = FALSE
  cohort$variantIDs.driver = unique(cohort$dataset[cohort$dataset$is.driver, 'variantID'])

  if(!is.null(cohort$phylogenies))
  {
    pio::pioTit("Propagating modification across trees stored in the cohort")

    # remove the driver even from the local copy of the data
    for(p in names(cohort$phylogenies))
    {
      for(f in 1:length(cohort$phylogenies[[p]]))
      {
        cohort$phylogenies[[p]][[f]]$dataset[
          cohort$phylogenies[[p]][[f]]$dataset$variantID %in% toDelete, 'is.driver'
        ] = FALSE

        # re-compute the information transfer
        cohort$phylogenies[[p]][[f]]$transfer = information.transfer(cohort$phylogenies[[p]][[f]])
      }
    }

  }

  # cat(cyan('\nChecking cohort.\n'))
  revolver_check_cohort(cohort)

  return(cohort)
}





#' @title Plot data and trees for a REVOLVER cohort.
#'
#' @details
#' Iterative plotting functions that scans a cohort and runs \code{\link{revolver_report_patient}}
#' on a set of patients (default all).
#'
#' @param x An object of class \code{"rev_cohort"}.
#' @param patients The patients to plot, default is all the one available.
#' @param max.phylogenies How many trees should be computed for each patient.
#' @param cex Scale cex for graphics.
#'
#' @return nothing
#' @export plot.rev_cohort
#' @import crayon
#'
#' @examples
#' data(Breast.fit)
#' plot.rev_cohort(Breast.fit, patients = Breast.fit$patients[1:10])
# plot.rev_cohort = function(x,
#                            patients = x$patients,
#                            max.phylogenies = 12,
#                            cex = 1)
# {
#   obj_has_trees(x)
#   plot.stat = TRUE
#
#   pio::pioHdr('REVOLVER Plot: Cohort (models)',
#               c(
#                 `Patients`=paste(patients, collapse = ', '),
#                 `Number of trees per patient` = max.phylogenies),
#               prefix = '\t -')
#
#   if(is.na(file)) stop('A file is required for this plot!')
#
#   for (patient in patients)
#   {
#     pio::pioTit(paste("Processing", patient))
#
#     revolver_report_patient(x, patient, cex = cex, max.phylogenies = max.phylogenies)
#   }
# }


