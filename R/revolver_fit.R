revolver_fit = function(x,
                        initial.solution = 1,
                        max.iterations = 10,
                        restarts = 0,
                        parallel = FALSE,
                        cores.ratio = .8,
                        transitive.orderings = FALSE,
                        verbose = FALSE)
{

  pio::pioHdr(paste0("REVOLVER fit - ", x$annotation))

  # What we can actually fit
  numPatients = length(x$phylogenies)
  fitPatients = names(x$phylogenies)

  # =-=-=-=-=-=-=-=-=-=-=-
  # Check if the cohort can be fit, stop on error
  # =-=-=-=-=-=-=-=-=-=-=-

  # Check types etc.
  revolver_check_cohort(x, stopOnError = TRUE)

  # Check input patients
  if(numPatients <= 1)
    stop("Cannot fit a model unless there are multiple patients with available trees, aborting")

  # Check restarts
  if(restarts > 0 & !is.na(initial.solution))
  {
    message('Beware: because you have set a fixed initial condition `restarts` will be disregarded because this EM is exhaustive.')

    restarts = 0
    # paralllel = FALSE
  }

  # =-=-=-=-=-=-=-=-=-=-=-
  # Print some info on the fitting task
  # =-=-=-=-=-=-=-=-=-=-=-

  # Fitting target
  pio::pioStr('\nFitting ', paste0('n = ', numPatients, 'patients'), suffix = '\n\n')
  print(Stats_trees(x, fitPatients))

  # Initial condition
  pio::pioStr('\nInitial solution :', ifelse(!is.na(initial.solution), paste0('Fixed to #', initial.solution), 'Randomized (uniform probability)'), suffix = '\n')
  pio::pioStr('\nSampled solutions (restarts) :', restarts, suffix = '\n')

  will_run_parallel = getOption("easypar.parallel", default = NA)
  pio::pioStr('\nParallel execuion (via \'easypar\') :',
              ifelse(is.na(will_run_parallel), TRUE, will_run_parallel)
              , suffix = '\n')

  results = NULL

  # Parallel fit
  if(parallel) {

    pclusters = setup_parallel(cores.ratio)

    r = foreach(num = 0:restarts, .packages = c("crayon", "igraph", 'matrixStats', 'matrixcalc'), .export = ls(globalenv())) %dopar%
    {
      tl_revolver_fit(x, initial.solution, max.iterations, transitive.orderings, verbose)
    }

    results = r

    stop_parallel(pclusters)
  }
  else {
    # NON Parallel fit
    # Fancy prints
    pline = function() {cat('\n\n\t'); for(i in 1:20)cat((cyan('~/-\\~')))}

    for(i in 0:restarts) {
      pline()
      model_fit = tl_revolver_fit(x, initial.solution, max.iterations, transitive.orderings, verbose)
      results = append(results, list(model_fit))
    }
  }

  # Get the best solution
  if(restarts > 0)
  {
    pline()
    cat(cyan('\n\nMedian goodness-of-fit (from fit penalty): '))
    scores = sapply(results, function(w) median(w$fit$penalty))
    cat(paste(scores, collapse = ', '))
    cat(cyan('\t Best is'), bgGreen('',scores[which.max(scores)], ''), '\n')
    return(results[[which.max(scores)]])
  }
  else return(results[[1]])
}



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# These are the functions that make the actual TL fit with an EM
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


############################## TL main function
#' @importFrom stats sd
tl_revolver_fit = function(x,
                           initial.solution = 1,
                           max.iterations = 10,
                           transitive.orderings = TRUE,
                           verbose = FALSE)
{
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Scoring function for a tree, in vectorized form
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  # Returns the score of a tree for a multinomial penalty E, and
  # and alpha parameter. It can also return just the penalty (i.e.,
  # without multiplying it for the actual score of the tree structure)
  Penalize = function(x, patient, E, rank = 1, alpha = 1, just.pen = FALSE) {
    IT = ITransfer(x, patient, rank)

    prod_penalties = IT %>%
      left_join(E, by = c("from", "to")) %>%
      pull(penalty) %>%
      prod()

    if(is.na(prod_penalties)) prod_penalties = 0

    if(just.pen) return(prod_penalties^alpha)

    Phylo(x, patient, rank)$score * (prod_penalties^alpha)
  }
  Penalize = Vectorize(Penalize, vectorize.args = 'rank')

  ################################################ Some functions which make easier to implement TL
  # add.to.Matrix = function(M, df){
  #   for(i in 1:nrow(df))
  #     M[df[i, 'from'], df[i, 'to']] = M[df[i, 'from'], df[i, 'to']] + 1
  #   M
  # }
  #
  # del.from.Matrix = function(M, df){
  #   for(i in 1:nrow(df))
  #     M[df[i, 'from'], df[i, 'to']] = M[df[i, 'from'], df[i, 'to']] - 1
  #   M
  # }
  #
  # asWeightedMatrix = function(df){
  #   entries.names = unique(unlist(df[ c(1,2)]))
  #
  #   M = matrix(0, ncol = length(entries.names), nrow = length(entries.names))
  #   colnames(M) = rownames(M) = entries.names
  #
  #   for(i in 1:nrow(df))
  #     M[df[i, 'from'], df[i, 'to']] = M[df[i, 'from'], df[i, 'to']] + 1
  #   M
  # }
  #
  # substitute.group.with.graph = function(g, G, S, verbose = FALSE)
  # {
  #   g = as.character(g)
  #
  #   S.w = S
  #   S[S>1] = 1
  #
  #   if(verbose) {
  #     cat(yellow('* Substitute '), g, '\n')
  #     cat(cyan('* Adj Matrix\n'))
  #     print(MatrixToDataFrame(G))
  #     cat(cyan('* Substitution Matrix\n'))
  #     print(S)
  #     # print(MatrixToDataFrame(S))
  #     # stop('sss')
  #   }
  #
  #   parent = pi(G, g)
  #   childrens = children(G, g)
  #
  #   # print(parent)
  #   # print(childrens)
  #   # print(length(childrens))
  #
  #   # loops are bastards, we need a special modified root and leaves function
  #   spec.root = function(w){
  #     r.G = igraph::graph_from_adjacency_matrix(w)
  #     r.w = root(w)
  #
  #     if(!igraph::is_dag(r.G)) {
  #       if(verbose) cat(red('Loops detected -- really bad temporal orderings within a cluster ... processing connected components ... \n', paste(MatrixToEdges(w), collapse = ' ')))
  #       comp = names(igraph::components(r.G)$membership)
  #       r.w = c(r.w, comp)
  #       r.w = unique(r.w)
  #     }
  #     r.w
  #   }
  #
  #   spec.leaves = function(w){
  #     r.G = igraph::graph_from_adjacency_matrix(w)
  #     r.w = leaves(w)
  #
  #     if(!igraph::is_dag(r.G)) {
  #       if(verbose) cat(red('Loops detected -- really bad temporal orderings within a cluster ... processing connected components ... \n', paste(MatrixToEdges(w), collapse = ' ')))
  #       comp = names(igraph::components(r.G)$membership)
  #       r.w = c(r.w, comp)
  #       r.w = unique(r.w)
  #     }
  #     r.w
  #   }
  #
  #   roots = spec.root(S)
  #   leaves = spec.leaves(S)
  #
  #   # print(roots)
  #   # print(leaves)
  #   attach.up = attach.down = NULL
  #   attach.up = apply(
  #     expand.grid(parent, roots, stringsAsFactors = FALSE),
  #     1,
  #     function(w)
  #       paste(w[1], w[2], sep = '~'))
  #   if(length(childrens) > 0) attach.down = apply(
  #     expand.grid(leaves, childrens, stringsAsFactors = FALSE),
  #     1,
  #     function(w)
  #       paste(w[1], w[2], sep = '~'))
  #
  #   # print(S)
  #   # print(roots)
  #   # print(expand.grid(parent, roots, stringsAsFactors = FALSE))
  #
  #   # paste(leaves, childrens, sep = '~')
  #
  #   delete.up = delete.down = NULL
  #   delete.up = paste(parent, g, sep = '~')
  #   if(length(childrens) > 0) delete.down = paste(g, childrens, sep = '~')
  #
  #   if(verbose) {
  #     cat(cyan('* Delete edges (inward)'), delete.up, '\n')
  #     cat(cyan('* Delete edges (outward)'), delete.down, '\n')
  #     cat(cyan('* Attachment edges (inward)'), attach.up, '\n')
  #     cat(cyan('* Attachment edges (outward)'), attach.down, '\n')
  #   }
  #
  #   edges = c(attach.up, attach.down)
  #   if(!all(S == 0)) edges = c(edges, MatrixToEdges(S))
  #
  #   # print(edges)
  #   # print(S)
  #
  #   G = MatrixToEdges(G)
  #
  #   # print(G)
  #
  #   #
  #   # print(attach.up)
  #   # print(attach.down)
  #
  #   G = setdiff(G, c(delete.up, delete.down))
  #   G = c(G, edges)
  #
  #   return(edgesToMatrix(G))
  # }
  #
  # wrapTS = function(M){
  #   tryCatch(
  #     {
  #       TS = igraph::topo_sort(igraph::graph_from_adjacency_matrix(M), mode = 'out')$name
  #       return(TS)
  #     },
  #     warning = function(w) { },
  #     error = function(w) { },
  #     finally = {return(TS)}
  #   )
  # }

  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # What we fit
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  numPatients = length(x$phylogenies)
  fitPatients = names(x$phylogenies)

  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Initial condition
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  if(!is.na(initial.solution)) {
    solutionID = rep(initial.solution, numPatients)
    names(solutionID) = fitPatients
  }
  else{
    solutionID = unlist(lapply(x$phylogenies, function(p) sample(1:length(p), 1)))
    names(solutionID) = fitPatients
  }

  # bestScores = stopCondition = penalty =  rep(NA, numPatients)
  # names(bestScores) = names(stopCondition) = names(penalty) = fitPatients

  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Drivers that will be correlated, add GL
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  all.drivers = lapply(fitPatients, Drivers, x = x)
  all.drivers = Reduce(bind_rows, all.drivers) %>% pull(variantID) %>% unique()
  all.drivers = c(all.drivers, "GL")

  # df.penalty = data.frame(matrix(nrow = 0, ncol = 3))
  # colnames(df.penalty) = c('Step', 'VariantID', 'penalty')

  # tb_penalty = data.frame(
  #   VariantID = all.drivers,
  #   Step = 0,
  #   Penalty = NA,
  #   stringsAsFactors = FALSE
  # ) %>% as_tibble()

  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # A data structure that we keep to store the fit progress
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  tb_solutionID = Stats_trees(x, fitPatients)
  tb_solutionID$Solution = solutionID

  #
  #   data.frame(
  #   patientID = fitPatients,
  #   Step = 0,
  #   Penalty = NA,
  #   stringsAsFactors = FALSE
  # ) %>% as_tibble()

  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Step 1) Fit via EM of the best tree for each patient
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  pio::pioStr('\n\n1] Expectation Maximization ', '\n')

  cat(inverse(sprintf("\n%27s", "Number of Solutions")))
  sapply(tb_solutionID$numTrees, function(p) cat(paste(sprintf('%4s', p, '|', sep =''))))

  cat(inverse(sprintf("\n%27s", "Combinations of Transfer")))
  sapply(tb_solutionID$combInfTransf, function(p) cat(paste(sprintf('%4s', p, '|', sep =''))))

  cat(inverse(sprintf("\n\n%27s", "Initialization")))
  sapply(tb_solutionID$Solution, function(s) cat(paste(green(sprintf('%4s', s)), '|', sep ='')))
  cat('\n')

  numIter = 1
  E = NULL
  repeat{

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    # E-step
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    cat(bgBlue("\n#", sprintf('%-3d', numIter), ' : '))
    cat(cyan("\tE: "))

    E = E_step(x, tb_solutionID)
    cat(green('OK'), cyan("\tM: "))

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    # M-step
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    tb_solutionID$newSolution = M_step(x, E, st_trees)
    tb_solutionID = tb_solutionID %>%
      mutate(
        converged = Solution == newSolution,
        Solution = newSolution
        ) %>%
      select(-newSolution)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    # Print to screen some nice colouring etc.
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    for(i in 1:nrow(tb_solutionID))
    {
      string = sprintf('%4s |', tb_solutionID$Solution[i])

      if(tb_solutionID$converged[i]) cat(green(string))
      else cat(red(string))
    }

    # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    # Stopping conditions (convergence or interrupt)
    # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    if(all(tb_solutionID$converged)) break;
    if(numIter == max.iterations)
    {
      cat(red('\n\n == Interrupted -- reached number of max.iterations requested. =='))
      break
    }

    numIter = numIter + 1
  }

  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Report the penalty to screen
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  solutionID = tb_solutionID %>% pull(Solution)

  tb_solutionID$penalty = sapply(
    seq(solutionID),
    function(w) Penalize(x, tb_solutionID$patientID[w],
                         E, rank = solutionID[w],
                         alpha = 1, just.pen = TRUE))

  print(tb_solutionID)



  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Store fits into a new object
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  fit = list()

  fit$penalty = E
  fit$phylogenies = list()

  for (patient in fitPatients)
    fit$phylogenies[[patient]] = Phylo(x,
                                       patient,
                                       tb_solutionID %>% filter(id == !!patient) %>% pull(Solution))

  stop("qui")

  #######################################################################
  ################################################ TL FIT -- step 2

  # Transitive closure of all information transfers
  cat(cyan('\n\n2] Transfering orderings across patients : \n\n'))

  cat((sprintf('   %20s', 'Transitivities are')), ifelse(transitive.orderings, green('enabled'), red('disabled')), '\n\n')

  cat(cyan(sprintf('   %20s : ', 'Information Transfer')))
  fit$tcInfTransf = lapply(fit$phylogenies, information.transfer, transitive.closure = transitive.orderings, indistinguishable = FALSE)
  names(fit$tcInfTransf) = names(fit$phylogenies)

  fit$tcOrderings = Reduce(rbind, lapply(fit$tcInfTransf, function(w) w$drivers))
  fit$tcOrderings = asWeightedMatrix(fit$tcOrderings)
  fit$solutionID = solutionID
  cat(green('OK\n'))

  # Assign what computed so far to 'x'
  x$fit = fit

  # verbose = TRUE

  # Then we take all ML orderings, for every patient
  cat(cyan(sprintf('   %20s : ', 'Expansions')), '\n')
  x$fit$substitutions = lapply(fitPatients, transfer.orderings, x = x, verbose = verbose)
  names(x$fit$substitutions) = fitPatients

  # print(x$fit$substitutions)

  npatients = length(x$fit$phylogenies)
  pb = txtProgressBar(min = 0, max = npatients, style = 3)
  pb.status = getOption('revolver.progressBar', default = TRUE)


  x$fit$exploded = list()
  for(patient in 1:npatients)
  {
    # update progress bar
    if(pb.status) setTxtProgressBar(pb, patient)

    patient = names(x$fit$phylogenies)[patient]

    phylo = x$fit$phylogenies[[patient]]
    subst = x$fit$substitutions[[patient]]

    # get groups, and the adj matrix
    groups = unique(phylo$dataset[phylo$dataset$is.driver, 'cluster'])
    adj_matrix = phylo$adj_mat

    # visit groups in order according to a topological sort, ensures correct expansions
    TS = wrapTS(adj_matrix)
    TS = TS[TS %in% groups]

    for(g in TS) adj_matrix = substitute.group.with.graph(g, G = adj_matrix, S = subst[[g]], verbose = verbose)

    x$fit$exploded[[patient]] = adj_matrix
  }


  # Now that we have expanded, we have to update the information transfer
  x$fit$transfer = lapply(names(x$fit$phylogenies),
                          function(w){
                            information.transfer_exploded_nodes(x, patient = w, transitive.orderings)
                          } )
  names(x$fit$transfer) = names(x$fit$phylogenies)


  class(x) = "rev_cohort_fit"

  return(x)
}

E_step = function(x, tb_solutionID)
{
  # Obtain all information transfers from the current solutions
  E = apply(
    data.frame(tb_solutionID),
    1,
    function(entry)
    {
      ITransfer(x, p = entry['patientID'], rank = as.numeric(entry['Solution']), type = 'drivers')
    }
  )

  E = Reduce(bind_rows, E)

  # Count and normalize them
  E %>%
    group_by(from, to) %>%
    summarise(count = n()) %>%
    group_by(to) %>%
    mutate(penalty = count / sum(count)) %>%
    ungroup()
}

M_step = function(x, E, st_trees)
{
  apply(st_trees,
        1,
        function(w) {
          scores = Penalize(x, w['patientID'], E, rank = 1:w['numTrees'])
          which.max(scores)
        })
}


# Penalize(x, patient, E, rank = 1:3)
