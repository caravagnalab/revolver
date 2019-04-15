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
  compute_trees_denovo = any(is.null(precomputed.trees))

  pio::pioHdr(paste("REVOLVER Construct phylogenetic trees for", patient),
              c(
                `Use precomputed trees` = compute_trees_denovo,
                `Maximum state space size to switch from exhaustive to Montecarlo search` = options$sspace.cutoff,
                `Number of Montecarlo samples, if not exhaustive` = options$n.sampling,
                `Overwrite the tree if it is already available` = options$overwrite,
                `Maximumum number of trees to store, if multiple are available` = options$store.max
              ),
              prefix = '\t'
  )

  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  # This function will not run for certain combination of parameters
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  # Prepare output
  if(!has_patient_trees(x)) x$phylogenies = NULL

  # Check if you do not want to overwrite already computed trees
  if(has_patient_trees(x, patient) & !as.logical(options['overwrite'])) {

    message('Trees already available for', patient, 'and `options$overwrite = FALSE`.
            The input cohort will be returned without modifications.\n')

    return(x)
  }

  # At this point some sort of computaiton will be done, and trees and their
  # scores will be stored in these two variables (lists)
  TREES = SCORES = NULL

  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # This function allows to input precomputed trees and in that case
  # no actual computation is done by REVOLVER and the trees are just
  # used as if they were computed by the tool
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  if(!compute_trees_denovo)
  {
    pio::pioTit('Precomputed trees given as input')

    are_suitable_precomputed_trees(x, patient, precomputed.trees, precomputed.scores)

    TREES = precomputed.trees
    SCORES = precomputed.scores
  }

  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # This branch computes de novo all the trees for a patient
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  if(compute_trees_denovo)
  {
    # We will need to use the CCF clusters for this patient
    clusters = CCF_clusters(x, patient)
    nclusters = nrow(clusters)

    clonal.cluster = Clonal_cluster(x, patient)
    samples = Samples(x, patient)

    df_clusters = data.frame(
      clusters[, samples],
      row.names = clusters$cluster
      )

    pio::pioTit(paste("CCF clusters for patient", patient))
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
      pio::pioTit("Building phylogentic trees, per region (this might take some time)")

      # clonal.cluster = as.character(unique(x$dataset$cluster[x$dataset$is.clonal]))
      # Get the data that we need to build a phylogenetic tree for a patient

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

      pio::pioDisp(data.frame(region = names(clonevol.obj$models), Trees = numSol, stringsAsFactors = FALSE))

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

      pio::pioTit("Scoring and ranking trees")

      # ################## Ranking trees. A tree is good according to the following factors:
      # # 1) the MI among the variables x and y, if they are connected by an edge x --> y [TODO: consider if we really need MI]
      # # 2) the Multinomial probability of edge x --> y in the trees determined by the CCF
      # # 3) for every edge  x --> y, the number of times that the CCF of x is greater than the CCF of y
      # # 3) for every node  x --> y1 ... yK, the number of times that the CCF of x is greater than the sum of the CCFs of y1 ... yK
      # binary.data = revolver:::binarize(x$dataset, samples)

      # 1) MI from binarized data -- different options, with a control sample which avoids 0log(0)
      # • a=0:maximum likelihood estimator (see entropy.empirical)
      # • a=1/2:Jeffreys’ prior; Krichevsky-Trovimov (1991) entropy estimator
      # • a=1:Laplace’s prior
      # • a=1/length(y):Schurmann-Grassberger (1996) entropy estimator
      # • a=sqrt(sum(y))/length(y):minimax prior
      binary.data = t(df_clusters)
      binary.data[binary.data > 0] = 1

      MI.table = revolver:::computeMI.table(binary.data, MI.Bayesian.prior = 0, add.control = TRUE)
      # Old parameter now fixed to FALSE
      use.MI = FALSE
      if(!use.MI) MI.table[TRUE] = 1

      # Steps 1 and 2 are collapsed, multiply MI by the Multinomial probability
      MI.table = revolver:::weightMI.byMultinomial(MI.table, WEIGHTS.CONSENSUS.TREE)

      # 3) Get penalty for direction given CCFs -- this is done for all possible edges in the data
      # CCF = clusters[, samples, drop = FALSE]
      # penalty.CCF.direction = edge.penalty.for.direction(TREES, CCF)
      penalty.CCF.direction = 1

      # 4) Compute the branching penalty  --  this is done for each tree that we are considering
      cat('Computing Pigeonhole Principle\n')
      penalty.CCF.branching = revolver:::node.penalty.for.branching(TREES, df_clusters)

      cat('Computing rank\n')
      RANKED = revolver:::rankTrees(TREES, MI.table, penalty.CCF.branching)
      TREES = RANKED$TREES
      SCORES = RANKED$SCORES

      TREES = TREES[SCORES > 0]
      SCORES = SCORES[SCORES > 0]

      TREES = lapply(TREES, DataFrameToMatrix)
    }
  }

  # =-=-=-=-=-=-=-=-=-=-=-=-=-=
  # Now we can transform the tree in our own format
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=
  pio::pioTit("Creating `rev_phylo`` object for REVOLVER (this might take some time)")

  # Special case: there are none
  if(length(TREES) == 0)
  {
    message(
      'No trees for this patient -- check data and input? -- returning original cohort.\n'
      )

    return(x)
  }

  # Now we create them
  x$phylogenies[[patient]] = create_trees_in_revolver_format(
    options,
    TREES,
    SCORES,
    patient,
    x,
    samples)

  # Just show how many combinations we have
  comb = combination_of_information_transfer(x, patient)
  pio::pioStr('\n Combinations of Information Transfer : ', comb, suffix = '\n')

  return(x)
}


# Check format of precomputed trees
are_suitable_precomputed_trees = function(
  x,
  patient,
  precomputed.trees = NULL,
  precomputed.scores = NULL
)
{
  # N clusters
  required_clusters = CCF_clusters(x, patient)$cluster

  # Check if they're all N x N
  cmt_size = sapply(precomputed.trees, ncol)
  rmt_size = sapply(precomputed.trees, nrow)

  OK_size = (cmt_size == rmt_size) && (cmt_size == length(required_clusters))

  if(!all(OK_size)) {
    message("Input trees must be NxN adjacency matrices because this patient has N = ",
            length(required_clusters), " clusters. Some are not.")

    stop("Cannot use these trees, aborting")
  }

  # Check if they have all N clusters
  OK_clusters = sapply(
    precomputed.trees,
    function(w) {
      all(required_clusters %in% colnames(w)) &&
        all(required_clusters %in% rownames(w))
    }
  )

  if(!all(OK_clusters)) {
    message("Input trees must contain entries for all clusters ",
            paste(required_clusters, collapse = ', '), " but some are not.")

    stop("Cannot use these trees, aborting")
  }

  # Check if they're actual trees
  OK_trees = sapply(precomputed.trees, is_tree)

  if(!all(OK_trees)) {
    message("Input trees are not trees and might have multiple roots, disconnected components etc.")

    stop("Cannot use these trees, aborting")
  }

  # All scores are not NA
  OK_scores_na = !(sapply(precomputed.scores, is.na))

  if(!all(OK_scores_na)) {
    message("Some trees are scored as NA, while they should all be real values.")

    stop("Cannot use these trees, aborting")
  }

  # Sorted scores
  OK_scores_sort = !(is.unsorted(rev(precomputed.scores)))

  if(!all(OK_scores_sort)) {
    message("Trees should be passed in descreasing order of scores, but they are not.")

    stop("Cannot use these trees, aborting")
  }
}
