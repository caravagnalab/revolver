#' @title Compute/ add CCF-based phylogenies to a REVOLVER cohort.
#'
#' @details
#'
#' Create or compute phylogenies for REVOLVER fit from CCF data. This method should be used also if one
#' has already pre-computed trees to use for fitting.
#'
#' @param x An object of class \code{"rev_cohort"}
#' @param patient A patient in the cohort, for which phylogenies are created
#' @param precomputed.trees If a list a of precomputed trees is available, this list should contain their
#' adjacency matrix. No computation will be carried out in this case.
#' @param precomputed.scores If a list a of precomputed trees is available, this list should contain their
#' scores. No computation will be carried out in this case.
#' @param options If one wants to generate new phylogenies for this patient, there named parameters should be
#' passed through this list. \code{sspace.cutoff = 10000} is the state-space cutoff to generate trees in an
#' exhaustive fashion, or to Montecarlo sample them. If Montecarlo is chosen \code{sampling} tells how
#' many trees are sampled. \code{overwrite = FALSE} sets if the patient's current phylogenies
#' should be overwritten, in case they are already available. \code{store.max = 100} tells how many trees should
#' be stored, if there are more than \code{store.max} available (ranked).
#' @param verbose output type.
#'
#' @return a modififed object of class \code{"rev_cohort"} with available phylogeneis for \code{patient}.
#' @export
#' @import crayon
#'
#' @examples
#' \dontrun{
#'  TODO
#' }
#'
revolver_compute_phylogenies = function(
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
  # Old parameter now fixed to FALSE
  use.MI = FALSE

  # Prepare output
  if(is.null(x$phylogenies)) x$phylogenies = NULL


  pio::pioHdr(paste("REVOLVER Construct phylogenetic trees for", patient),
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

  pCCF = CCF(x, patient)
  samples = Samples(x, patient)

  dataset = left_join(Data(x, patient), pCCF, by = 'id')

  # cat(cyan('[compute_rev_phylogenies] Clusters for this patient ... \n'))
  clusters = revolver:::clusters.table(x$dataset, samples)
  nclusters = nrow(clusters)

  TREES = SCORES = NULL

  if(!any(is.null(precomputed.trees))){
    cat(yellow('\n\nPrecomputed trees given as input ... using them.\n\n'))

    TREES = precomputed.trees
    SCORES = precomputed.scores
  }
  else
  {
    pio::pioTit("Groups/ Clusters in the data of this patient")
    print(clusters)

    if(nclusters == 1)
    {
      cat(red('\nThis model has 1 node, it has trivial models.'))

      M = matrix(0, ncol = 1, nrow = 1)
      colnames(M) = rownames(M) = rownames(clusters)

      TREES = append(TREES, list(M))
      SCORES = c(SCORES, 1)
    }
    else
    {
      # ################## Generate all trees that are compatible with the observed CCFs, we do this
      # ################## by analyzing one sample at a time.
      pio::pioTit("Using ClonEvol to build phylogentic trees, per region (this might take some time)")

      clonal.cluster = as.character(unique(x$dataset$cluster[x$dataset$is.clonal]))

      # December 2018 - https://github.com/caravagn/revolver/issues/13
      #
      # Original code which would run ClonEvol's functions
      #
      # x$dataset[, samples] =  x$dataset[, samples] * 100
      # clonevol.obj = revolver:::useClonevo(x$dataset, samples, clonal.cluster)
      # x$dataset[, samples] =  x$dataset[, samples] / 100
      #
      # Alternative implementation of ClonEvol's steps that we use. We use our function
      # and reacreate an obj matching the expected output, so we do not need
      # to update downstream functions.
      alternative = NULL
      alternative$models = ClonEvol_surrogate(clusters, samples, clonal.cluster, min.CCF = 0.01)
      clonevol.obj = alternative

      # remove the trees which have no edges (returned for samples with only 1 cluster for instance)
      numSol = sapply(clonevol.obj$models, function(w){ sum(sapply(w, nrow) > 0) })
      pio::pioDisp(data.frame(region = names(clonevol.obj$models), Solutions = numSol))

      ################## Build all possible clonal trees
      # 1) hash them
      # 2) create a consensus as the union of all trees
      # 3) generate or sample a large number of possible trees, where a parent x --> y is assigned
      #    with probability proportional to how often the edge is detected

      pio::pioTit("Hashing, merging and generating solutions (this might take some time)")
      CLONAL.TREES = revolver:::hashTrees(clonevol.obj, samples)

      CONSENSUS = revolver:::consensusModel(clonevol.obj, samples)
      CONSENSUS.TREE = CONSENSUS$S
      WEIGHTS.CONSENSUS.TREE = CONSENSUS$weights

      cat("Created graph of all possible ancestries\n")

      # # Sampling is carried out if there are more than 'sspace.cutoff' trees, in that case we
      # # sample 'n.sampling' possible trees. Otherwise all possible trees are generated.
      TREES = revolver:::all.possible.trees(
        CONSENSUS.TREE,
        WEIGHTS.CONSENSUS.TREE,
        sspace.cutoff = options['sspace.cutoff'],
        n.sampling = options['n.sampling']
      )

      pio::pioTit("Scoring and ranking solutions (this might take some time)")

      # ################## Ranking trees. A tree is good according to the following factors:
      # # 1) the MI among the variables x and y, if they are connected by an edge x --> y [TODO: consider if we really need MI]
      # # 2) the Multinomial probability of edge x --> y in the trees determined by the CCF
      # # 3) for every edge  x --> y, the number of times that the CCF of x is greater than the CCF of y
      # # 3) for every node  x --> y1 ... yK, the number of times that the CCF of x is greater than the sum of the CCFs of y1 ... yK
      binary.data = revolver:::binarize(x$dataset, samples)

      # 1) MI from binarized data -- different options, with a control sample which avoids 0log(0)
      # • a=0:maximum likelihood estimator (see entropy.empirical)
      # • a=1/2:Jeffreys’ prior; Krichevsky-Trovimov (1991) entropy estimator
      # • a=1:Laplace’s prior
      # • a=1/length(y):Schurmann-Grassberger (1996) entropy estimator
      # • a=sqrt(sum(y))/length(y):minimax prior
      MI.table = revolver:::computeMI.table(binary.data, MI.Bayesian.prior = 0, add.control = TRUE)
      if(!use.MI) MI.table[TRUE] = 1

      # Steps 1 and 2 are collapsed, multiply MI by the Multinomial probability
      MI.table = revolver:::weightMI.byMultinomial(MI.table, WEIGHTS.CONSENSUS.TREE)

      # 3) Get penalty for direction given CCFs -- this is done for all possible edges in the data
      CCF = clusters[, samples, drop = FALSE]
      # penalty.CCF.direction = edge.penalty.for.direction(TREES, CCF)
      penalty.CCF.direction = 1

      # 4) Compute the branching penalty  --  this is done for each tree that we are considering
      cat('Computing Pigeonhole Principle\n')
      penalty.CCF.branching = revolver:::node.penalty.for.branching(TREES, CCF)

      cat('Computing rank\n')
      RANKED = revolver:::rankTrees(TREES, MI.table, penalty.CCF.branching)
      TREES = RANKED$TREES
      SCORES = RANKED$SCORES

      TREES = TREES[SCORES > 0]
      SCORES = SCORES[SCORES > 0]
    }
  }

  if(length(TREES) == 0){
    cat(red('No phylogenies found for this patient -- check data? -- returning original cohort.\n'))

    x$dataset = Original.dataset
    return(x)
  }

  pio::pioTit("Creating rev_phylo object for REVOLVER (this might take some time)")

  x$phylogenies[[patient]] = create_trees_in_revolver_format(options, TREES, SCORES, patient, x$dataset, samples)

  # Restore data
  x$dataset = Original.dataset

  comb = rev_count_information_transfer_comb(x, patient)
  pio::pioStr('\n Combinations of Information Transfer : ', crayon::yellow(comb), suffix = '\n')

  return(x)
}
