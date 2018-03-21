#' Construct a REVOLVER cohort object (S3 class \code{"rev_cohort"}).
#'
#' @param dataset A dataframe in the specified format (see Online manual).
#' @param CCF.parser A function to parse the format for the encoding of CCF
#' or binary values for each sequenced region. A possible function is available
#' inside REVOLVER; since it is not exported but is available with
#' \code{revolver:::CCF.parser} (the default of this parameter).
#' @param options A list of 2 parameters that should be a boolean value for
#' \code{ONLY.DRIVER} (use only driver SNVs), and \code{MIN.CLUSTER.SIZE}, the minimum cluster size.
#' @param annotation A string for annotation of this cohort. This will be prompted
#' in every print for this object.
#'
#' @return An object of class \code{"rev_cohort"}
#' @export
#' @import crayon
#'
#' @examples
#' data(CRC)
#' cohort = revolver_cohort(CRC)
revolver_cohort = function(
  dataset,
  CCF.parser = revolver:::CCF.parser,
  options = list(ONLY.DRIVER = FALSE, MIN.CLUSTER.SIZE = 10),
  annotation = '')
{
  REVOLVER_VERSION = "\"Haggis and tatties\""
  # by G.Caravagna <giulio.caravagna@icr.ac.uk>"

  # Dataset with certain columns
  required.cols = c('Misc', 'patientID', 'variantID',  'cluster', 'is.driver', 'is.clonal', 'CCF')
  if(!all(required.cols %in% colnames(dataset)))
    stop('Dataset should have the following columns:', paste(required.cols, collapse = ', '))

  if(!is.function(CCF.parser))
     stop('You need to provide a function to parse CCFs.')

  cat(bgCyan('REVOLVER'), cyan(REVOLVER_VERSION))
  cat(blue(paste(' -- ', annotation, sep ='')), '\n')

  cat('Options:\n')
  cat('\t ONLY.DRIVER:', options$ONLY.DRIVER, '\n')
  cat('\t MIN.CLUSTER.SIZE:', options$MIN.CLUSTER.SIZE, '\n')

  dataset = dataset[, required.cols]

  if(options$ONLY.DRIVER)
  {
    cat('\n[DRIVER = F] Dropping all non DRIVER variants ... ')
    dataset = dataset[dataset$is.driver, , drop = FALSE]
    cat('OK\n')
  }

  if(options$MIN.CLUSTER.SIZE > 0){
    cat('\n[MIN.CLUSTER.SIZE > 0] filtering clusters below that size ... \n')

    cat('\tNumber of rows before filtering:', nrow(dataset), '\n')

    data.split = split(dataset, f = dataset$patientID)

    data.reduced = lapply(
      data.split,
      function(p){
        p = filter.clusters(p, cutoff.numMuts = options$MIN.CLUSTER.SIZ)
      })

    dataset = Reduce(rbind, data.reduced)
    cat('\tNumber of rows after filtering:', nrow(dataset), '\n')
  }

  patients = unique(dataset$patientID)
  variantIDs = unique(dataset$variantID)
  variantIDs.driver = unique(dataset[dataset$is.driver, 'variantID'])
  numVariants = nrow(dataset)

  obj <-
    structure(
      list(
        patients = patients,
        variantIDs = variantIDs,
        variantIDs.driver = variantIDs.driver,
        numVariants = numVariants,
        annotation = annotation,
        dataset = dataset,
        CCF.parser = CCF.parser,
        REVOLVER_VERSION = REVOLVER_VERSION
      ),
      class = "rev_cohort",
      call = match.call()
    )

  return(obj)
}

#' Print a \code{"rev_cohort"} object
#'
#' @param x obj of class \code{"rev_cohort"}
#' @param digits number of output digits
#'
#' @return nothing
#' @export
#' @import crayon
#'
#' @examples
#' data(CRC.cohort)
#' CRC.cohort
print.rev_cohort <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  stopifnot(inherits(x, "rev_cohort"))

  cat(bgCyan('REVOLVER'), cyan(x$REVOLVER_VERSION))
  cat(blue(paste(' -- ', x$annotation, sep ='')), '\n')

  cat(cyan('\n\tPatients      :'), paste(head(x$patients), collape =', ', sep =''), '... \n')
  cat(cyan('\tNum. Patients :'), length(x$patients), '\n')
  cat(cyan('\tVariants      :'), x$numVariants, 'of which', length(x$variantIDs.driver), 'drivers\n')

  cat(cyan('\n\tPhylogenies   :'),
      ifelse(is.null(x$phylogenies), red('NO'), green('YES')))

  cat(cyan('\n\tFit           :'),
      ifelse(is.null(x$fit), red('NO'), green('YES')))

  cat(cyan('\n\tClusters      :'),
      ifelse(is.null(x$cluster), red('NO'), green('YES')))


  if(!is.null(x$phylogenies))
  {
    # cat('-- There are phylogenies for', length(names(x$phylogenies)), 'patient(s)')

    cat(cyan('\n\tTL Model Fit  :'), ifelse(is.null(x$fit), red('NO\n'), green('YES ')), '\n')

    longest.name = max(nchar(names(x$phylogenies)))

    patf = function(w, fit){
      s = paste(
        yellow(sprintf(paste('%', longest.name, 's', sep = ''), names(x$phylogenies)[w])),  ': ',
          sprintf('k = %3s | t = %3s | n = %2s | r = %2s | m = %4s | d = %2s',
                  length(x$phylogenies[[w]]),
                  rev_count_information_transfer_comb(x, names(x$phylogenies)[w]),
                  x$phylogenies[[w]][[1]]$numNodes,
                  x$phylogenies[[w]][[1]]$numRegions,
                  nrow(x$phylogenies[[w]][[1]]$dataset),
                  nrow(x$phylogenies[[w]][[1]]$dataset[x$phylogenies[[w]][[1]]$dataset$is.driver, ])
          ),  '\t')

      cat('\t', s)

      if(fit)
      {
        cat(bgBlue('[ Fit ]'))
        stas = stats.rev_phylo(x$fit$phylogenies[[w]])
        cat(sprintf(' # %3s', blue(x$fit$solutionID[w])))
        cat('| g', red(sprintf('%8s', round(stas$gofit, 6))))
        cat('| f', green(sprintf('%s', x$fit$phylogenies[[w]]$score)))
      }
      cat('\n')
    }

    sapply(1:length(x$phylogenies), patf, fit = !is.null(x$fit))

    cat('\n\tLegend \n\t\t k : phylogenies  \n\t\t t : combinations of information transfer \n\t\t n : groups (nodes of the tree) \n\t\t r : regions (inputs per patient) \n\t\t m : number of alterations \n\t\t d : number of driver alterations\n')

    if(!is.null(x$fit))
      cat('\t\t # : number of the solution selection (out of k)  \n\t\t g : goodness-of-fit \n\t\t f : score of the model\n')

  }

  cat('\n')
  revolver_check_cohort(x)

}

rev_count_information_transfer_comb = function(x, p) {
  stopifnot(!is.null(x$phylogenies))
  stopifnot(p %in% names(x$phylogenies))

  keys = lapply(
    x$phylogenies[[p]],
        function(w)
          paste(sort(DataFrameToEdges(w$transfer$driver)), collapse = ' '))

  keys = Reduce(rbind, keys)
  length(unique(keys))
}



#' @title Compute/ add CCF-based phylogenies to a REVOLVER cohort.
#'
#' @details
#'
#' Create or compute phylogenies for REVOLVER fit from CCF data. This method should be used also if one
#' has already pre-computed trees to use for fitting.
#'
#' @param x An object of class \code{"rev_cohort"}
#' @param patient A patient in the cohort, for which phylogenies are created
#' @param use.MI If MI should be used to weight each phylogenetic model (should be FALSE).
#' @param precomputed.trees If a list a of precomputed trees is available, this list should contain their
#' adjacency matrix. No computation will be carried out in this case.
#' @param precomputed.trees If a list a of precomputed trees is available, this list should contain their
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
#' @examples TODO
revolver_compute_phylogenies = function(
  x,
  patient,
  use.MI = FALSE,
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

  cat(cyan('\n\nCreating phylogenies for'), yellow(patient))

  if(!is.null(x$phylogenies))
  {
    if(patient %in% names(x$phylogenies))
    {
      if(as.logical(options['overwrite']))
        cat(red('\t overwrite = TRUE -- Overwriting them.'))
      else
        {
          cat(red('overwrite = FALSE -- Skipping patient.'))
          Sys.sleep(0.5)
          return(x)
        }
    }
  }
  cat('\n')



  # Phylo creation
  x$dataset = x$dataset[x$dataset$patientID  == patient, ]
  # cat(cyan('[compute_rev_phylogenies] Entries for this patient: '), nrow(x$dataset), '\n')

  samples = names(x$CCF.parser(x$dataset[1, 'CCF']))

  # CCF = Reduce(rbind, sapply(x$dataset$CCF, x$CCF.parser))
  # cat(cyan('CCF '))

  CCF = sapply(x$dataset$CCF, x$CCF.parser)
  # print(CCF)
  if(length(samples) == 1) {
    CCF = matrix(CCF, ncol = 1)
    colnames(CCF) = samples
  } else CCF = t(CCF)

  CCF = apply(CCF, 2, as.numeric)
  rownames(CCF) = rownames(x$dataset)

  # cat(green('OK'))

  # cat(cyan('[compute_rev_phylogenies] Binding CCF -- samples:'), paste(colnames(CCF), collapse = ', '), '\n')
  x$dataset = x$dataset[, c('Misc', 'patientID', 'variantID', 'cluster', 'is.driver', 'is.clonal')]
  x$dataset = cbind(x$dataset, CCF)

  if(verbose) print(head(x$dataset))

  # cat(cyan('[compute_rev_phylogenies] Clusters for this patient ... \n'))
  clusters = clusters.table(x$dataset, samples)
  nclusters = nrow(clusters)

  print(clusters)


  TREES = SCORES = NULL
  if(!any(is.null(precomputed.trees))){
    cat(yellow('\nPrecomputed trees given as input ... using them.\n'))

    TREES = precomputed.trees
    SCORES = precomputed.scores
  }
  else
  {
    # cat(yellow('\nPhylogenies for', patient, ':'))

    if(nclusters == 1)
    {
      cat(red('1 single cluster (no-edge model with fictitious score 1)'))

      M = matrix(0, ncol = 1, nrow = 1)
      colnames(M) = rownames(M) = rownames(clusters)

      TREES = append(TREES, list(M))
      SCORES = c(SCORES, 1)
    }
    else
    {
      # ################## Generate all trees that are compatible with the observed CCFs, we do this
      # ################## by analyzing one sample at a time.
      # cat(yellow('1) Generate all trees that are compatible with the observed CCFs ...\n'))

      x$dataset[, samples] =  x$dataset[, samples] * 100

      cat(cyan('Clonevo: '))
      clonal.cluster = as.character(unique(x$dataset$cluster[x$dataset$is.clonal]))
      clonevol.obj = useClonevo(x$dataset, samples, clonal.cluster)

      x$dataset[, samples] =  x$dataset[, samples] / 100


      cat(paste(names(clonevol.obj$models), unlist(lapply(clonevol.obj$models, length)), collapse = ', '), ' ')


      ################## Build all possible clonal trees
      # 1) hash them
      # 2) create a consensus as the union of all trees
      # 3) generate or sample a large number of possible trees, where a parent x --> y is assigned
      #    with probability proportional to how often the edge is detected
      # cat(yellow('2) Build or sample all possible clonal trees ...\n'))

      # cat(cyan('[compute_rev_phylogenies] Hashing trees\n'))
      cat(cyan(' | Trees '))
      CLONAL.TREES = hashTrees(clonevol.obj, samples)
      # print(CLONAL.TREES)

      # cat(cyan('[compute_rev_phylogenies] Computing consensus\n'))
      cat(cyan('| Consensus '))
      CONSENSUS = consensusModel(clonevol.obj, samples)
      CONSENSUS.TREE = CONSENSUS$S
      WEIGHTS.CONSENSUS.TREE = CONSENSUS$weights
      cat('OK')


      # # Sampling is carried out if there are more than 'sspace.cutoff' trees, in that case we
      # # sample 'n.sampling' possible trees. Otherwise all possible trees are generated.
      # cat(cyan('[compute_rev_phylogenies] Generating possible solutions (it might take time) ... \n'))
      cat(cyan('\nSolutions '))
      TREES = all.possible.trees(
        CONSENSUS.TREE,
        WEIGHTS.CONSENSUS.TREE,
        sspace.cutoff = options['sspace.cutoff'],
        n.sampling = options['n.sampling']
      )

      # ################## Ranking trees. A tree is good according to the following factors:
      # # 1) the MI among the variables x and y, if they are connected by an edge x --> y [TODO: consider if we really need MI]
      # # 2) the Multinomial probability of edge x --> y in the trees determined by the CCF
      # # 3) for every edge  x --> y, the number of times that the CCF of x is greater than the CCF of y
      # # 3) for every node  x --> y1 ... yK, the number of times that the CCF of x is greater than the sum of the CCFs of y1 ... yK
      #       ################## Generate binary data from CCFs
      binary.data = binarize(x$dataset, samples)
      # cat(cyan('[compute_rev_phylogenies] Binary data generated\n'))

      # cat(yellow('3) Ranking trees.\n'))

      # 1) MI from binarized data -- different options, with a control sample which avoids 0log(0)
      # • a=0:maximum likelihood estimator (see entropy.empirical)
      # • a=1/2:Jeffreys’ prior; Krichevsky-Trovimov (1991) entropy estimator
      # • a=1:Laplace’s prior
      # • a=1/length(y):Schurmann-Grassberger (1996) entropy estimator
      # • a=sqrt(sum(y))/length(y):minimax prior
      MI.table = computeMI.table(binary.data, MI.Bayesian.prior = 0, add.control = TRUE)
      if(!use.MI) MI.table[TRUE] = 1

      # print(MI.table)
      # print(WEIGHTS.CONSENSUS.TREE)

      # Steps 1 and 2 are collapsed, multiply MI by the Multinomial probability
      MI.table = weightMI.byMultinomial(MI.table, WEIGHTS.CONSENSUS.TREE)

      # cat(cyan('[compute_rev_phylogenies] Mutual Information weighted by Multinomial score\n'))
      # print(MI.table)
      cat(cyan('Mutual Information'), ifelse(use.MI, green('YES'), red('NO')))



      # 3) Get penalty for direction given CCFs -- this is done for all possible edges in the data
      CCF = clusters[, samples, drop = FALSE]
      # penalty.CCF.direction = edge.penalty.for.direction(TREES, CCF)
      penalty.CCF.direction = 1

      # 4) Compute the branching penalty  --  this is done for each tree that we are considering
      # cat(cyan('[compute_rev_phylogenies] Computing violations of the Pigeonhole principle for every model\n'))
      cat(cyan(' | Pigeonhole principle \n'))
      penalty.CCF.branching = node.penalty.for.branching(TREES, CCF)



      # cat(cyan('[compute_rev_phylogenies] Ranking trees, and removing 0-scoring ones\n'))
      cat(cyan('Ranking trees\n'))
      RANKED = rankTrees(TREES, MI.table, penalty.CCF.branching)
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


  if(length(TREES) > options['store.max']){
    cat(cyan('Phylogenies found'), length(TREES), red('-- storing', options['store.max'], ' '))
    TREES = TREES[1:as.numeric(options['store.max'])]
    SCORES = SCORES[1:as.numeric(options['store.max'])]
  }

  cat(cyan('Creating'), length(TREES), cyan('revolver_phylogeny objects\n'))

  LSCORES = as.data.frame(SCORES)
  LSCORES = split(LSCORES, f = LSCORES[,1])
  # print(SCORES)
  LSCORES = lapply(LSCORES, function(w) w[sample(1:nrow(w), replace = FALSE), , drop = FALSE])
  # print(SCORES)

  permuted.indexes = as.integer(rev(unlist(lapply(LSCORES, rownames))))
  names(permuted.indexes) = NULL
  #
  #
  # print(SCORES)
  # print(LSCORES)
  # print(permuted.indexes)

  TREES = TREES[permuted.indexes]
  SCORES = SCORES[permuted.indexes]

  # print(TREES)
  # print(SCORES)

  # print(permuted.indexes)

  REVOLVER.TREES = NULL
  for(i in 1:length(TREES))
  {
    cat('@ ', i, '\r')

    tree = revolver_phylogeny(
      M = TREES[[i]],
      patient = patient,
      dataset = x$dataset,
      samples = samples,
      score = SCORES[i],
      annotation = paste('Ranked ', i, '/', length(TREES), sep ='')
    )

    # print('afaasas')

    REVOLVER.TREES = append(REVOLVER.TREES, list(tree))
  }
  x$dataset = Original.dataset
  if(is.null(x$phylogenies)) x$phylogenies = NULL
  x$phylogenies[[patient]] = REVOLVER.TREES

  comb = rev_count_information_transfer_comb(x, patient)
  cat(cyan('#Information Transfer'), ifelse(comb == 1, red(comb), green(comb)), '\n\n')



  return(x)

  # ################## Generate all trees that are compatible with the observed CCFs, we do this
  # ################## by analyzing one sample at a time.
  # clonal.cluster = as.character(unique(my.data$cluster[my.data$is.clonal]))
  # clonevol.obj = useClonevo(my.data, sample.groups, clonal.cluster)
  #
  # STAT.ENTRY = c(
  #   STAT.ENTRY,
  #   numPhylogenies = paste(names(clonevol.obj$models), unlist(lapply(clonevol.obj$models, length)), collapse = ', '))
  #
  # ################## Generate binary data from CCFs
  # binary.data = binarize(my.data, sample.groups)
  #
  # ################## Build all possible clonal trees
  # # 1) hash them
  # # 2) create a consensus as the union of all trees
  # # 3) generate or sample a large number of possible trees, where a parent x --> y is assigned
  # #    with probability proportional to how often the edge is detected
  # CLONAL.TREES = hashTrees(clonevol.obj, sample.groups)
  #
  # CONSENSUS = consensusModel(clonevol.obj, sample.groups)
  # CONSENSUS.TREE = CONSENSUS$S
  # WEIGHTS.CONSENSUS.TREE = CONSENSUS$weights
  #
  # # Sampling is carried out if there are more than 'sspace.cutoff' trees, in that case we
  # # sample 'n.sampling' possible trees. Otherwise all possible trees are generated.
  # TREES = all.possible.trees(
  #   CONSENSUS.TREE,
  #   WEIGHTS.CONSENSUS.TREE,
  #   sspace.cutoff = 10000,
  #   n.sampling = 5000
  # )
  #
  # ################## Ranking trees. A tree is good according to the following factors:
  # # 1) the MI among the variables x and y, if they are connected by an edge x --> y [TODO: consider if we really need MI]
  # # 2) the Multinomial probability of edge x --> y in the trees determined by the CCF
  # # 3) for every edge  x --> y, the number of times that the CCF of x is greater than the CCF of y
  # # 3) for every node  x --> y1 ... yK, the number of times that the CCF of x is greater than the sum of the CCFs of y1 ... yK
  #
  # # 1) MI from binarized data -- different options, with a control sample which avoids 0log(0)
  # # • a=0:maximum likelihood estimator (see entropy.empirical)
  # # • a=1/2:Jeffreys’ prior; Krichevsky-Trovimov (1991) entropy estimator
  # # • a=1:Laplace’s prior
  # # • a=1/length(y):Schurmann-Grassberger (1996) entropy estimator
  # # • a=sqrt(sum(y))/length(y):minimax prior
  # MI.table = computeMI.table(binary.data, MI.Bayesian.prior = 0, add.control = TRUE)
  #
  # # Steps 1 and 2 are collapsed, multiply MI by the Multinomial probability
  # MI.table = weightMI.byMultinomial(MI.table, WEIGHTS.CONSENSUS.TREE)
  #
  # # 3) Get penalty for direction given CCFs -- this is done for all possible edges in the data
  # CCF = clusters.table(my.data, sample.groups)[, sample.groups]
  # penalty.CCF.direction = edge.penalty.for.direction(TREES, CCF)
  #
  # # 4) Compute the branching penalty  --  this is done for each tree that we are considering
  # penalty.CCF.branching = node.penalty.for.branching(TREES, CCF)
  #
  # structural.score = apply(
  #   cbind(penalty.CCF.direction,penalty.CCF.branching),
  #   1,
  #   function(x) prod(x))
  #
  # RANKED = rankTrees(TREES, MI.table, structural.score)
  # TREES = RANKED$TREES
  # SCORES = RANKED$SCORES
  #
  # TREES = TREES[SCORES > 0]
  # SCORES = SCORES[SCORES > 0]
  #
  # STAT.ENTRY = c(
  #   STAT.ENTRY,
  #   numSolutions = length(TREES))
  #
  # REVOLVER.TREES = NULL
  # for(x in 1:length(TREES))
  # {
  #   cat('@ ', x, '\r')
  #
  #   tree = revolver_phylogeny(
  #     M = TREES[[x]],
  #     patient = PATIENT,
  #     dataset = my.data,
  #     score = SCORES[x],
  #     annotation = paste('Ranked ', x, '/', length(TREES), sep ='')
  #   )
  #
  #   REVOLVER.TREES = append(REVOLVER.TREES, list(tree))
  # }
  #
  # # quartz()
  # # lapply(REVOLVER.TREES, plot, plot.stat = T, file = 's.pdf')
  # #
  # plot(REVOLVER.TREES[[1]], edge.label = round(MI.table, 2))
  #
  # # TREES = weighted.sampling(
  # #   DataFrameToMatrix(CONSENSUS.TREE),
  # #   WEIGHTS.CONSENSUS.TREE,
  # #   1000
  # # )
  # # TREES = checkColinearity(TREES, binary.data)
  #
  # quartz()
  # plotCCFTrees(
  #   CLONAL.TREES,
  #   my.data,
  #   binary.data,
  #   CONSENSUS,
  #   MI.table,
  #   REVOLVER.TREES,
  #   SCORES,
  #   PATIENT,
  #   max.trees.in.plots = 10,
  #   cex = 1)
  #
  # ccf = list(
  #   CLONAL.TREES = CLONAL.TREES,
  #   my.data = my.data,
  #   binary.data = binary.data,
  #   CONSENSUS = CONSENSUS,
  #   MI.table = MI.table,
  #   TREES = TREES,
  #   PATIENT = PATIENT
  # )
  #
  # save(ccf,file = paste(PATIENT, '-CCF-Analysis.RData'))
  #
  # PROCESSED = c(PROCESSED, PATIENT)
  # STATS = rbind(STATS, STAT.ENTRY)
  #
  # cat(PATIENT, 'processed ... moving to next patient in 2 seconds.\n')
  # Sys.sleep(2)

}

#' @title  Compute/ add CCF-based mutation trees to a REVOLVER cohort.
#'
#' @details
#' This is the analogous of \code{\link{revolver_compute_phylogenies}}, but for binary data and
#' hence it computes Chow-Liu trees using also Suppes' conditions. Parameters have exactly
#' the same meaning of the ones described in \code{\link{revolver_compute_phylogenies}}.
#'
#' @param x An object of class \code{"rev_cohort"}
#' @param patient A patient in the cohort, for which mutation trees are created
#' @param precomputed.trees See \code{\link{revolver_compute_phylogenies}}
#' @param precomputed.trees See \code{\link{revolver_compute_phylogenies}}
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

  cat(cyan('\n\nCreating phylogenies for'), yellow(patient))

  if(!is.null(x$phylogenies))
  {
    if(patient %in% names(x$phylogenies))
    {
      if(as.logical(options['overwrite']))
        cat(red('\t overwrite = TRUE -- Overwriting them.'))
      else
      {
        cat(red('overwrite = FALSE -- Skipping patient.'))
        Sys.sleep(0.5)
        return(x)
      }
    }
  }
  cat('\n')



  # Phylo creation
  x$dataset = x$dataset[x$dataset$patientID  == patient, ]
  # cat(cyan('[compute_rev_phylogenies] Entries for this patient: '), nrow(x$dataset), '\n')

  samples = names(x$CCF.parser(x$dataset[1, 'CCF']))

  # CCF = Reduce(rbind, sapply(x$dataset$CCF, x$CCF.parser))
  # cat(cyan('CCF '))

  CCF = sapply(x$dataset$CCF, x$CCF.parser)
  # print(CCF)
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


  # cat(cyan('[compute_rev_phylogenies] Binding CCF -- samples:'), paste(colnames(CCF), collapse = ', '), '\n')
  x$dataset = x$dataset[, c('Misc', 'patientID', 'variantID', 'cluster', 'is.driver', 'is.clonal')]
  x$dataset = cbind(x$dataset, CCF)

  if(verbose) print(head(x$dataset))

  # cat(cyan('[compute_rev_phylogenies] Clusters for this patient ... \n'))
  clusters = clusters.table(x$dataset, samples)
  nclusters = nrow(clusters)

  print(clusters)


  TREES = SCORES = NULL
  if(!any(is.null(precomputed.trees))){
    cat(yellow('\nPrecomputed trees given as input ... using them.\n'))

    TREES = precomputed.trees
    SCORES = precomputed.scores
  }
  else
  {
    # cat(yellow('\nPhylogenies for', patient, ':'))

    if(nclusters == 1)
    {
      cat(red('1 single cluster (no-edge model with fictitious score 1)'))

      M = matrix(0, ncol = 1, nrow = 1)
      colnames(M) = rownames(M) = rownames(clusters)

      TREES = append(TREES, list(M))
      SCORES = c(SCORES, 1)
    }
    else
    {
      ################### Generate Suppes poset
      cat(cyan('| Suppes '))
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
      cat('OK')

      # # Sampling is carried out if there are more than 'sspace.cutoff' trees, in that case we
      # # sample 'n.sampling' possible trees. Otherwise all possible trees are generated.
      # cat(cyan('[compute_rev_phylogenies] Generating possible solutions (it might take time) ... \n'))
      cat(cyan('\nSolutions '))
      TREES = all.possible.trees(
        G = POSET.EDGES,
        W = POSET,
        sspace.cutoff = options['sspace.cutoff'],
        n.sampling = options['n.sampling']
      )

      # ################## Ranking trees. A tree is good according to the following factors:
      # # 1) the MI among the variables x and y, if they are connected by an edge x --> y [TODO: consider if we really need MI]
      # # 2) the Multinomial probability of edge x --> y in the trees determined by the CCF
      # # 3) for every edge  x --> y, the number of times that the CCF of x is greater than the CCF of y
      # # 3) for every node  x --> y1 ... yK, the number of times that the CCF of x is greater than the sum of the CCFs of y1 ... yK
      #       ################## Generate binary data from CCFs
      binary.data = binarize(x$dataset, samples)
      # cat(cyan('[compute_rev_phylogenies] Binary data generated\n'))

      # cat(yellow('3) Ranking trees.\n'))

      # 1) MI from binarized data -- different options, with a control sample which avoids 0log(0)
      # • a=0:maximum likelihood estimator (see entropy.empirical)
      # • a=1/2:Jeffreys’ prior; Krichevsky-Trovimov (1991) entropy estimator
      # • a=1:Laplace’s prior
      # • a=1/length(y):Schurmann-Grassberger (1996) entropy estimator
      # • a=sqrt(sum(y))/length(y):minimax prior
      MI.table = computeMI.table(binary.data, MI.Bayesian.prior = 0, add.control = TRUE)
      cat(cyan('Mutual Information'), green('YES '))

      cat(cyan('Ranking trees\n'))
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


  if(length(TREES) > options['store.max']){
    cat(cyan('Phylogenies found'), length(TREES), red('-- storing', options['store.max'], ' '))
    TREES = TREES[1:as.numeric(options['store.max'])]
    SCORES = SCORES[1:as.numeric(options['store.max'])]
  }

  cat(cyan('Creating'), length(TREES), cyan('revolver_phylogeny objects\n'))

  LSCORES = as.data.frame(SCORES)
  LSCORES = split(LSCORES, f = LSCORES[,1])
  # print(SCORES)
  LSCORES = lapply(LSCORES, function(w) w[sample(1:nrow(w), replace = FALSE), , drop = FALSE])
  # print(SCORES)

  permuted.indexes = as.integer(rev(unlist(lapply(LSCORES, rownames))))
  names(permuted.indexes) = NULL
  #
  #
  # print(SCORES)
  # print(LSCORES)
  # print(permuted.indexes)

  TREES = TREES[permuted.indexes]
  SCORES = SCORES[permuted.indexes]

  # print(TREES)
  # print(SCORES)

  # print(permuted.indexes)

  REVOLVER.TREES = NULL
  for(i in 1:length(TREES))
  {
    cat('@ ', i, '\r')

    tree = revolver_phylogeny(
      M = TREES[[i]],
      patient = patient,
      dataset = x$dataset,
      samples = samples,
      score = SCORES[i],
      annotation = paste('Ranked ', i, '/', length(TREES), sep ='')
    )

    # print('afaasas')

    REVOLVER.TREES = append(REVOLVER.TREES, list(tree))
  }
  x$dataset = Original.dataset
  if(is.null(x$phylogenies)) x$phylogenies = NULL
  x$phylogenies[[patient]] = REVOLVER.TREES

  comb = rev_count_information_transfer_comb(x, patient)
  cat(cyan('#Information Transfer'), ifelse(comb == 1, red(comb), green(comb)), '\n\n')



  return(x)

  # ################## Generate all trees that are compatible with the observed CCFs, we do this
  # ################## by analyzing one sample at a time.
  # clonal.cluster = as.character(unique(my.data$cluster[my.data$is.clonal]))
  # clonevol.obj = useClonevo(my.data, sample.groups, clonal.cluster)
  #
  # STAT.ENTRY = c(
  #   STAT.ENTRY,
  #   numPhylogenies = paste(names(clonevol.obj$models), unlist(lapply(clonevol.obj$models, length)), collapse = ', '))
  #
  # ################## Generate binary data from CCFs
  # binary.data = binarize(my.data, sample.groups)
  #
  # ################## Build all possible clonal trees
  # # 1) hash them
  # # 2) create a consensus as the union of all trees
  # # 3) generate or sample a large number of possible trees, where a parent x --> y is assigned
  # #    with probability proportional to how often the edge is detected
  # CLONAL.TREES = hashTrees(clonevol.obj, sample.groups)
  #
  # CONSENSUS = consensusModel(clonevol.obj, sample.groups)
  # CONSENSUS.TREE = CONSENSUS$S
  # WEIGHTS.CONSENSUS.TREE = CONSENSUS$weights
  #
  # # Sampling is carried out if there are more than 'sspace.cutoff' trees, in that case we
  # # sample 'n.sampling' possible trees. Otherwise all possible trees are generated.
  # TREES = all.possible.trees(
  #   CONSENSUS.TREE,
  #   WEIGHTS.CONSENSUS.TREE,
  #   sspace.cutoff = 10000,
  #   n.sampling = 5000
  # )
  #
  # ################## Ranking trees. A tree is good according to the following factors:
  # # 1) the MI among the variables x and y, if they are connected by an edge x --> y [TODO: consider if we really need MI]
  # # 2) the Multinomial probability of edge x --> y in the trees determined by the CCF
  # # 3) for every edge  x --> y, the number of times that the CCF of x is greater than the CCF of y
  # # 3) for every node  x --> y1 ... yK, the number of times that the CCF of x is greater than the sum of the CCFs of y1 ... yK
  #
  # # 1) MI from binarized data -- different options, with a control sample which avoids 0log(0)
  # # • a=0:maximum likelihood estimator (see entropy.empirical)
  # # • a=1/2:Jeffreys’ prior; Krichevsky-Trovimov (1991) entropy estimator
  # # • a=1:Laplace’s prior
  # # • a=1/length(y):Schurmann-Grassberger (1996) entropy estimator
  # # • a=sqrt(sum(y))/length(y):minimax prior
  # MI.table = computeMI.table(binary.data, MI.Bayesian.prior = 0, add.control = TRUE)
  #
  # # Steps 1 and 2 are collapsed, multiply MI by the Multinomial probability
  # MI.table = weightMI.byMultinomial(MI.table, WEIGHTS.CONSENSUS.TREE)
  #
  # # 3) Get penalty for direction given CCFs -- this is done for all possible edges in the data
  # CCF = clusters.table(my.data, sample.groups)[, sample.groups]
  # penalty.CCF.direction = edge.penalty.for.direction(TREES, CCF)
  #
  # # 4) Compute the branching penalty  --  this is done for each tree that we are considering
  # penalty.CCF.branching = node.penalty.for.branching(TREES, CCF)
  #
  # structural.score = apply(
  #   cbind(penalty.CCF.direction,penalty.CCF.branching),
  #   1,
  #   function(x) prod(x))
  #
  # RANKED = rankTrees(TREES, MI.table, structural.score)
  # TREES = RANKED$TREES
  # SCORES = RANKED$SCORES
  #
  # TREES = TREES[SCORES > 0]
  # SCORES = SCORES[SCORES > 0]
  #
  # STAT.ENTRY = c(
  #   STAT.ENTRY,
  #   numSolutions = length(TREES))
  #
  # REVOLVER.TREES = NULL
  # for(x in 1:length(TREES))
  # {
  #   cat('@ ', x, '\r')
  #
  #   tree = revolver_phylogeny(
  #     M = TREES[[x]],
  #     patient = PATIENT,
  #     dataset = my.data,
  #     score = SCORES[x],
  #     annotation = paste('Ranked ', x, '/', length(TREES), sep ='')
  #   )
  #
  #   REVOLVER.TREES = append(REVOLVER.TREES, list(tree))
  # }
  #
  # # quartz()
  # # lapply(REVOLVER.TREES, plot, plot.stat = T, file = 's.pdf')
  # #
  # plot(REVOLVER.TREES[[1]], edge.label = round(MI.table, 2))
  #
  # # TREES = weighted.sampling(
  # #   DataFrameToMatrix(CONSENSUS.TREE),
  # #   WEIGHTS.CONSENSUS.TREE,
  # #   1000
  # # )
  # # TREES = checkColinearity(TREES, binary.data)
  #
  # quartz()
  # plotCCFTrees(
  #   CLONAL.TREES,
  #   my.data,
  #   binary.data,
  #   CONSENSUS,
  #   MI.table,
  #   REVOLVER.TREES,
  #   SCORES,
  #   PATIENT,
  #   max.trees.in.plots = 10,
  #   cex = 1)
  #
  # ccf = list(
  #   CLONAL.TREES = CLONAL.TREES,
  #   my.data = my.data,
  #   binary.data = binary.data,
  #   CONSENSUS = CONSENSUS,
  #   MI.table = MI.table,
  #   TREES = TREES,
  #   PATIENT = PATIENT
  # )
  #
  # save(ccf,file = paste(PATIENT, '-CCF-Analysis.RData'))
  #
  # PROCESSED = c(PROCESSED, PATIENT)
  # STATS = rbind(STATS, STAT.ENTRY)
  #
  # cat(PATIENT, 'processed ... moving to next patient in 2 seconds.\n')
  # Sys.sleep(2)

}


#' @title Plot a REVOLVER cohort.
#'
#' @details
#' Iterative plotting functions that scans samples from a cohort. It will plot
#' a patient's data (original and binarized), the score for the trees associated to
#' the patient, and a number of possible trees per patient.
#'
#' @param x An object of class \code{"rev_cohort"}
#' @param patients The patients to plot, default is all the one available.
#' @param max.phylogenies How many trees should be computed for each patient
#' @param cex Scale cex for graphics
#' @param plot.stat Set it to TRUE to visualize also some statistics about each tree
#'
#' @return nothing
#' @export
#' @import RColorBrewer
#' @import pheatmap
#' @import TRONCO
#' @import crayon
#'
#' @examples
#' data(CRC.cohort)
#' plot(CRC.cohort, patients = 'adenoma_3')
plot.rev_cohort = function(x, patients = x$patients, max.phylogenies = 12, cex = 1, plot.stat = TRUE)
{
  for(patient in patients)
  {
    cat(yellow(patient), ':')

    ################################# Clusters table
    cat('\tCCF')

    pat = x$phylogenies[[patient]][[1]]
    ccf = pat$CCF[, pat$samples_ID, drop = F]

    ccf.numbers = ccf
    ccf.numbers = round(ccf.numbers, 3)
    ccf.numbers[ccf.numbers == 0] = ''

    annot = pat$CCF[, c('is.driver', 'is.clonal')]
    annot$is.clonal = as.character(annot$is.clonal)
    annot$is.driver = as.character(annot$is.driver)

    # annot = apply(annot, 2, as.character)

    yes_no_col = c('forestgreen', 'gainsboro')
    names(yes_no_col) = c('TRUE', 'FALSE')
    ann_colors = list(is.driver = yes_no_col, is.clonal = yes_no_col)

    # pheatmap(ccf)

    br = c(0, 1e-3, seq(0.1, 1, 0.1))
    capture.output({pheatmap::pheatmap(
      ccf,
      main = 'CCF',
      cluster_cols = FALSE,
      cluster_rows = FALSE,
      border = NA,
      breaks = br,
      legend = F,
      color = c('white', colorRampPalette(brewer.pal(9, 'Blues'))(length(br) - 1)),
      display_numbers = ccf.numbers,
      cellwidth = 34,
      cellheight = 12,
      number_color = 'orange',
      # width = 3,
      annotation_row = annot,
      annotation_colors = ann_colors,
      # height = 8, width = 8,
      # height = 8,
      file = paste(patient, '-CCF.pdf', sep = '')
    )})

    ################################# Binarized data
    cat('| Binary data')

    capture.output({
    TRONCO::oncoprint(
      TRONCO::import.genotypes(rbind(pat$binary.data, wt = 0), color = 'steelblue'),
      cellwidth = 22,
      cellheight = 10,
      sample.id = TRUE,
      font.row = 8,
      ann.hits = FALSE,
      legend = FALSE,
      title = 'Binarized data',
      file = paste(patient, '-data.pdf', sep = '')
    )})

    ################################# Scores for phylogenies
    cat('| Phylogeny scores')

    scores = lapply(x$phylogenies[[patient]], function(e) e$score)
    comb = rev_count_information_transfer_comb(x, patient)
    colors = colorRampPalette(brewer.pal(8, 'Dark2'))(comb)

    groups = lapply(x$phylogenies[[patient]],
                    function(w) paste(sort(DataFrameToEdges(w$transfer$drivers)), collapse = ':'))
    groups = unlist(groups)
    names(colors) = unique(groups)

    pdf(paste(patient,'-scoreplot.pdf', sep =''), width = 2, height = 3)
    barplot(unlist(scores), col = colors[groups], border = NA, main = patient,  horiz = T, ylab = 'Phylogenies', xlab = 'Score')
    legend('topright', legend = comb, bty = 'n')
    dev.off()

    ##### Plot the distribution
    distribution = x$phylogenies[[patient]]

    if(length(distribution) > max.phylogenies) distribution = distribution[1:max.phylogenies]
    cat('| ', max.phylogenies, ' Phylogenies :')


    toMerge = sapply(
      1:length(distribution),
      function(x){
        f = paste(patient, x, 'distribution.pdf', sep = '-')

        plot(
          distribution[[x]],
          plot.stat = plot.stat,
          graph.cex = cex,
          table.cex = cex,
          # edge.label = round(MI.table, 2),
          file = f
        )

        cat('+')

        return(f)
      })

    xx = jamPDF(
      in.files = toMerge,
      out.file =  paste(patient, '-distribution.pdf', sep = ''),
      layout = '2x2'
    )

    xx = jamPDF(in.files = paste(patient, c('-CCF.pdf', '-data.pdf', '-scoreplot.pdf'), sep =''), out.file = paste(patient, '-scoreplot.pdf', sep = ''), layout = '3x1')

    xx =  jamPDF(
      in.files = paste(patient, c('-scoreplot.pdf', '-distribution.pdf'), sep =''),
      out.file =  paste(patient, '.pdf', sep = ''),
      layout = '1x1'
    )

    cat('\n')
  }
}


clonal.subclonal.table = function(x)
{
  clonal.drivers = table(x$dataset[x$dataset$is.driver & x$dataset$is.clonal, 'variantID'])
  subclonal.drivers = table(x$dataset[x$dataset$is.driver & !x$dataset$is.clonal, 'variantID'])

  tableD = data.frame(
    variantID = unique(x$dataset[x$dataset$is.driver, 'variantID']),
    stringsAsFactors = FALSE)
  rownames(tableD) = tableD$variantID

  tableD = cbind(tableD, Clonal = 0)
  tableD = cbind(tableD, SubClonal = 0)

  for(d in tableD$variantID){
    tableD[d, 'Clonal'] = clonal.drivers[d]
    tableD[d, 'SubClonal'] = subclonal.drivers[d]
  }
  tableD[is.na(tableD)] = 0
  tableD$Counts = tableD$Clonal + tableD$SubClonal

  tableD = tableD[order(tableD$Counts, decreasing = T), ]
  tableD$variantID = NULL

  tableD
}


#' @title Check basic inconsistencies in a REVOLVER cohort.
#'
#' @details Perform some basic diagnostic of a cohort object. It will inform of patients without drivers and other
#' information that can be used to reshape the data before fitting a model.
#'
#' @param x An object of class \code{"rev_cohort"}
#' @param auto.fix Unused parameter; ideally we will implement some automatic fixing one day.
#' @param return.value If the method should return or not a report with the errors. By default is just prints to screen
#'
#' @return depends on \code{return.value}
#' @import crayon
#' @export
#'
#' @examples
#' data(CRC.cohort)
#' revolver_check_cohort(CRC.cohort)
#' print(CRC.cohort) # calls this anyway
revolver_check_cohort = function(x, auto.fix = FALSE, return.value = FALSE)
{
  # check for duplicated IDs ... parallel evolution events
  data = split(x$dataset, f = x$dataset$patientID)

  err = FALSE
  xx = lapply(
    split(x$dataset, f = x$dataset$patientID),
    function(w)
    {
      df = w[w$is.driver, ]
      mask = duplicated(df$variantID)

      if(any(mask))
      {
        cat(bgRed('\n CHECK '), red("Duplicated driver variantID in patient: "), yellow(df$patientID[1]), ':', df$variantID[mask], '\n')
        entries = df[df$variantID %in% df$variantID[mask], ]

        err = TRUE

        print(entries)
        return(entries)
      }
    })

  # Check for variants cocurring <= 1 times
  data.split = x$dataset[x$dataset$is.driver, ]
  data.split = split(data.split, f = data.split$variantID)

  occurrencesCount = lapply(data.split, function(x) nrow(x) <= 1)
  if(any(unlist(occurrencesCount))) {
    cat(bgRed('\n CHECK '), red("These drivers occurr once, will not be correlated: "))
    occurrencesCount = occurrencesCount[unlist(occurrencesCount)]

    err = TRUE

    cat(paste(names(occurrencesCount), collapse = ', '))
    # for(o in names(occurrencesCount)) {
    #   cat(yellow(o), '@', data.split[[o]][1, 'patientID'], ' | ')
    #
    #   # if(auto.fix) {
    #   #   cat(green('AUTOFIX [= TRUE] will delete it.'))
    #   #   df = data.split[[o]]
    #   #
    #   #   x = revolver_remove.driver(x,  df['patientID'], df['variantID'], df['Misc'], df['cluster'])
    #   #
    #   # }
    # }
    # cat('\n')
  }

  # check for patients with no variants
  data.split = split(x$dataset, f = x$dataset$patientID)
  data.split = lapply(data.split, function(w) w[w$is.driver, ])

  zoesp = lapply(data.split, function(w) nrow(w) == 0)
  unoesp = lapply(data.split, function(w) nrow(w) == 1)

  if(any(unlist(zoesp))) {
    cat(bgRed('\n CHECK '), red("Patients with 0 drivers are useless and should be removed:"), paste(names(zoesp)[unlist(zoesp)], collapse = ', '))

    err = TRUE
  }

  if(any(unlist(unoesp)))
    cat(bgRed('\n CHECK '), red("Patients with 1 drivers can only be expanded:"), paste(names(unoesp)[unlist(unoesp)], collapse = ', '))


  # if(any(unlist(occurrencesCount))) stop('Errors should be fixed.')

  # return(x)

  # if(auto.fix)
  # {
  #   cat(green('\nAUTOFIX running on these entries.\n'))
  #   print(Reduce(rbind, xx))
  #
  #   for(i in 1:length(xx)){
  #   }
  # }

 if(return.value) return(err)
}

#' @title Remove a driver event from a cohort.
#'
#' @details Basic editing function. Each event is identied through its ID,
#' cluster assignment and its misc flag. With this function, you can remove it.
#' The event is not physically removed from the dataset, but instead its flag
#' \code{is.Driver} is set to \code{FALSE}. If there are phylogenies inside, they are updated
#' as well. If you have fit the models, however, you should re-run the fit after this
#' modification because this modification does not propagate.
#'
#' @param x An object of class \code{"rev_cohort"}
#' @param patient Patient ID.
#' @param variantID ID of the event to mark as \code{FALSE}.
#' @param misc misc of the event to mark as \code{FALSE}
#' @param cluster cluster of the event to mark as \code{FALSE}
#'
#' @return a modified cohort object of class \code{"rev_cohort"}
#' @import crayon
#' @export
#'
#' @examples
#' data(CRC.cohort)
#' head(CRC.cohort$dataset) # Info to use to remove
#' new.data = revolver_removeDriver(CRC.cohort, "adenoma_1", "APC", "NOTHING", "1")
#' head(new.data$dataset)
revolver_removeDriver = function(
  x,
  patient,
  variantID,
  misc,
  cluster)
{
  x$dataset[
    x$dataset$is.driver &
      x$dataset$patientID == patient &
      x$dataset$Misc == misc &
      x$dataset$cluster == cluster, 'is.driver'
  ] = FALSE

  if(!is.null(x$phylogenies)){
    for(p in names(x$phylogenies)){
      for(model in 1: length(x$phylogenies[[p]])) {
        x$phylogenies[[p]][[model]]$dataset[
          x$dataset$is.driver &
            x$dataset$patientID == patient &
            x$dataset$Misc == misc &
            x$dataset$cluster == cluster, 'is.driver'
          ] = FALSE
      }
    }

  }



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
  current.drivers = rownames(clonal.subclonal.table(cohort))
  toDelete = setdiff(current.drivers, list)

  cat(cyan('Removing:'), paste(toDelete, collapse = ', '))

  cohort$dataset[cohort$dataset$variantID %in% toDelete, 'is.driver'] = FALSE
  cohort$variantIDs.driver = unique(cohort$dataset[cohort$dataset$is.driver, 'variantID'])

  if(!is.null(cohort$phylogenies))
  {
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

  cat(cyan('\nChecking cohort.\n'))
  revolver_check_cohort(cohort)

  return(cohort)
}



#' Subset the patients in the cohort to match the input list.
#'
#' @param cohort An object of class \code{"rev_cohort"}
#' @param list A vector of patient IDs to subset the data to.
#'
#' @return A REVOLVER cohort with only patients in the \code{list}.
#' @export
#'
#' @examples
#' data(CRC.cohort)
#' CRC.cohort
#' revolver_deletePatients(CRC.cohort, "adenoma_2")
revolver_deletePatients = function(x, list)
{
  new.patients = setdiff(x$patients, list)
  x$dataset = x$dataset[x$dataset$patientID %in% new.patients, , drop = FALSE]
  x$patients = new.patients

  if(!is.null(x$phylogenies))
    x$phylogenies = x$phylogenies[intersect(new.patients, names(x$phylogenies))]

  cat(cyan('\nChecking cohort.\n'))
  revolver_check_cohort(x)


  return(x)
}

CCF.parser = function(x)
{
  tk = strsplit(x, ';')[[1]]
  tk = unlist(strsplit(tk, ':'))

  samples = tk[seq(1, length(tk), 2)]

  values = tk[seq(2, length(tk), 2)]
  names(values) = samples

  return(values)
}

