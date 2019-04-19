

revolver_fit = function(x,
                        initial.solution = 1,
                        max.iterations = 10,
                        n = 10,
                        ...)
{
  pio::pioHdr(paste0("REVOLVER Transfer Learning fit- ", x$annotation))
  stopifnot(n > 0)

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
  if(n > 1 & !is.na(initial.solution))
  {
    message('Beware: because you have set a fixed initial condition `n` will be disregarded because this EM is exhaustive.')

    n = 1
  }

  # =-=-=-=-=-=-=-=-=-=-=-
  # Print some info on the fitting task
  # =-=-=-=-=-=-=-=-=-=-=-

  # Fitting target
  pio::pioStr('\nFitting ', paste0('N = ', numPatients, ' patients'), suffix = '\n\n')
  print(Stats_trees(x, fitPatients))

  # Initial condition
  pio::pioStr(
    '\nInitial solution :',
    ifelse(
      !is.na(initial.solution),
      paste0('Fixed to #', initial.solution),
      'Randomized (uniform probability)'
    ),
    suffix = '\n'
  )
  pio::pioStr('\nSampled solutions: ', paste0('n = ', n), suffix = '\n')

  will_run_parallel = getOption("easypar.parallel", default = NA)
  pio::pioStr(
    '\nParallel exectuion (via \'easypar\') :',
    ifelse(is.na(will_run_parallel), TRUE, will_run_parallel),
    suffix = '\n'
  )

  results = easypar::run(
    FUN = function(w)
    {
      print(ls())

      tl_revolver_fit(x,
                      initial.solution = initial.solution,
                      max.iterations = max.iterations)

    },
    PARAMS = lapply(1:n, list),
    packages = c("crayon", "revolver", "tidygraph", "tidyverse", "igraph"),
    export = ls(globalenv(), all.names = TRUE)
  )

  if(easypar::numErrors(results))
  {
    message("Returned the following errors")
    lapply(results,
           function(w) {
             if(inherits(w, 'simpleError') | inherits(w, 'try-error'))
               print(w)
             })
  }

  # Get the best solution
  best = 1

  if(n > 1)
  {
    pio::pioTit("Selecting solution with minimal median penalty")

    # We take the max of the penalty variable (which is the same because at 1 the penalty is 0)
    scores = sapply(seq_along(results), function(w) {

      run_score = median(results[[w]]$fit$fit_table$penalty)

      pio::pioStr(paste0("Solution #", w), run_score, suffix = '\n')
      })

    best = which.max(scores)

    cat(cyan('\t Best solution is #'), bgGreen(best), '\n')
  }

  return(results[[best]])
}



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# These are the functions that make the actual TL fit with an EM
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


############################## TL main function
#' @importFrom stats sd
tl_revolver_fit = function(x,
                           initial.solution = 1,
                           max.iterations = 10
                           )
{
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Auxiliary functions for the fitting algorithm
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  # Scoring function for a tree, in vectorized form
  #
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

  # EM functions
  E_step = function(x, tb_solutionID, E = NULL)
  {
    # This function can also get a list in E input, in that case it will just
    # skip this step and straight away do the binding etc. This allows to use
    # this function in different phases of the algorithm
    if(all(is.null(E)))
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
    }

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

  # Return for every clone in a patient, a possible ordering of its events
  # based from the MLE estimate that we obtain from the multinomial model (counts)
  ML_clones_expansion = function(x, patient, rank, E)
  {
    # Get clones
    drivers = Phylo(x, patient, rank)$drivers %>%
      rename(to = variantID)
      # %>%
      # group_by(cluster) %>%
      # filter(n() > 1)

    # Algorithm to expand a single group
    expand_group = function(g)
    {
      g_drivers = drivers %>% filter(cluster == !!g)

      # the MLE estimator for the drivers in a group is obtained couting
      # only the orderings among these drivers. We re-normalize it as well
      ML_E = E %>%
        filter(
          from %in% g_drivers$to,
          to %in% g_drivers$to
        ) %>%
        group_by(to) %>%
        mutate(penalty = count / sum(count)) %>%
        ungroup()

      # Create a graph with the nodes of this expansion
      nodes_gp = data.frame(from = 'GL', to = g_drivers$to, stringsAsFactors = FALSE)

      # The nodes are called "cluster" in this tb_graph
      as_tbl_graph(bind_rows(nodes_gp, ML_E)) %>%
        filter(name != 'GL') %>%
        rename(cluster = name)
    }

    # Each group is expanded
    expansions = lapply(
      unique(drivers$cluster),
      expand_group
    )
    names(expansions) = unique(drivers$cluster)

    expansions
  }

  # For a given patient, it takes each one of the expanded clones after the second step of
  # transfer learning, and connects them to their descendant as of the information transfer.
  Attach_expanded = function(x, patient, rank, tb_graphs)
  {
    # Extract from a graph the nodes that are part of a loop
    # as defined from the Strongly Connected Components.
    loopy_nodes = function(M)
    {
      G = igraph::graph_from_adjacency_matrix(M)

      # Strongly Connected Components have loops by definition
      SCC = igraph::components(G, "strong")$membership

      as_tibble(data.frame(node = names(SCC), SCC = SCC, stringsAsFactors = FALSE)) %>%
        group_by(SCC) %>%
        filter(n() > 1) %>%
        pull(node)
    }

    # Get the clones to expand from the transfer
    clones = ITransfer(x, patient, rank, type = 'clones')

    # The topological sort of the transfer to understand in which order they shall be traversed
    clones_orderings = igraph::topo_sort(
      igraph::graph_from_adjacency_matrix(DataFrameToMatrix(clones)),
      mode = 'out'
      )$name

    # Get the drivers which provide the actual nodes to attach, and add a special map for GL
    # drivers = Phylo(x, patient, rank)$drivers
    # drivers = bind_rows(drivers,
    #                     data.frame(
    #                       variantID = 'GL',
    #                       cluster = "GL",
    #                       stringsAsFactors = FALSE
    #                     ))

    # The expanded graph will start from GL
    # current_expansion = clones %>%
    #   filter(from == "GL")

    # The initial condition is that nothing is expanded
    current_expansion = clones

    # Expand clones following the topological ordering
    for (clone in clones_orderings)
    {
      # GL is not expanded by definition, any other node is
      if (clone == 'GL')
        next

      # Get the adjaceny matrix of what we are expanding, and the list of edges as well
      clone_adjM = TidyGraphToMatrix(tb_graphs[[clone]])
      clone_df = TidyGraphToDataFrame(tb_graphs[[clone]])

      # First we compute the roots and leaves of this subgraph
      nodes_roots = root(clone_adjM)
      nodes_leaves = leaves(clone_adjM)

      # Get loopy nodes that are tricky to handle, we treat them as both roots/leaves
      lnodes = loopy_nodes(M = clone_adjM)
      nodes_roots = c(nodes_roots, lnodes)
      nodes_leaves = c(nodes_leaves, lnodes)

      # Then we take what are currently expanding, and the parent and children in there
      current_expansion_adjM = DataFrameToMatrix(current_expansion)
      parent_current_expansion = pi(current_expansion_adjM, clone)
      children_current_expansion = children(current_expansion_adjM, clone)

      # We now make the modification, we can essentially
      # 1) take the current expansion and remove all edges that involve the current clone
      # 2) attach all roots to the parental node of clone, in the graph that we are expanding,
      # 2) attach all leaves to the children node of clone, in the graph that we are expanding,
      # 2) copy all the rest (intermediate nodes) in the current expansion
      current_expansion = bind_rows(
        current_expansion %>%
          filter(from != clone, to != clone),
        rbind(
          expand.grid(
            from = parent_current_expansion,
            to = nodes_roots,
            stringsAsFactors = FALSE
          ),
          expand.grid(
            from = nodes_leaves,
            to = children_current_expansion,
            stringsAsFactors = FALSE
          ),
          clone_df
        ) %>%
          as_tibble()
      )
    }

    current_expansion
  }

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

  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Drivers that will be correlated, add GL
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  all.drivers = lapply(fitPatients, Drivers, x = x)
  all.drivers = Reduce(bind_rows, all.drivers) %>% pull(variantID) %>% unique()
  all.drivers = c(all.drivers, "GL")

  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # A data structure that we keep to store the fit progress in the first EM
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  tb_solutionID = Stats_trees(x, fitPatients)
  tb_solutionID$Solution = solutionID

  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Step 1) Fit via EM of the best tree for each patient
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  pio::pioTit('1] Expectation Maximization')

  cat(inverse(sprintf("\n%27s", "Number of Solutions")))
  sapply(tb_solutionID$numTrees, function(p) cat(paste(sprintf('%4s', p, '|', sep =''))))

  cat(inverse(sprintf("\n%27s", "Combinations of Transfer")))
  sapply(tb_solutionID$combInfTransf, function(p) cat(sprintf('%4s |', p)))

  cat(inverse(sprintf("\n\n%27s", "Initialization")))
  sapply(tb_solutionID$Solution, function(s) cat(green(sprintf('%4s |', s))))
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

    tb_solutionID$newSolution = M_step(x, E, tb_solutionID)
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
    function(w) Penalize(x,
                         tb_solutionID$patientID[w],
                         E,
                         rank = solutionID[w],
                         alpha = 1,
                         just.pen = TRUE))

  print(tb_solutionID)

  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Store fits into a new object which we later add to `x`
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  fit = list()

  fit$fit_table = tb_solutionID
  fit$penalty = E

  fit$phylogenies = list()

  for (patient in fitPatients)
    fit$phylogenies[[patient]] = Phylo(x,
                                       patient,
                                       tb_solutionID %>% filter(patientID == !!patient) %>% pull(Solution))

  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # 2) ML parent set transferred from the other patients
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  pio::pioTit('2] Transfering orderings across patients')

  fit$clones_expansions = fit$clones_to_expand = NULL

  # We track this with a progress bar
  pb = txtProgressBar(min = 0, max = length(fit$phylogenies), style = 3)
  pb.status = getOption('revolver.progressBar', default = TRUE)

  # For every patient we obtain the ML expansions for each patient, and then we
  # assemble that in an expanded graph. After this, we add the expanded graph
  # as the Information Transfer field of the driver events in each fit phylogeny,
  # and we change the annotation of the model to make it clear
  for (patient in seq_along(fitPatients))
  {
    # update progress bar
    if(pb.status) setTxtProgressBar(pb, patient)

    patient = fitPatients[patient]

    fit$clones_to_expand[[patient]] = ML_clones_expansion(
        x,
        patient,
        rank = tb_solutionID %>% filter(patientID == !!patient) %>% pull(Solution),
        E
      )

    fit$clones_expansions[[patient]] = Attach_expanded(
        x,
        patient,
        rank = tb_solutionID %>% filter(patientID == !!patient) %>% pull(Solution),
        tb_graphs = fit$clones_to_expand[[patient]]
      )

    fit$phylogenies[[patient]]$transfer$drivers =
      fit$clones_expansions[[patient]]

    fit$phylogenies[[patient]]$annotation =
      paste0(
        fit$phylogenies[[patient]]$annotation,
        " - Information Transfer expanded via Transfer Learning"
      )
  }
  names(fit$clones_expansions) = names(fit$clones_to_expand) = fitPatients

  # We can now get the penalty from the final fit output, re-using the E-step function
  # because we are still computing an expectation here
  fit$penalty = E_step(x, NULL, fit$clones_expansions)

  x$fit = fit

  # cat((sprintf('   %20s', 'Transitivities are')), ifelse(transitive.orderings, green('enabled'), red('disabled')), '\n\n')
#
#   cat(cyan(sprintf('   %20s : ', 'Information Transfer')))
#   fit$tcInfTransf = lapply(fit$phylogenies, information.transfer, transitive.closure = transitive.orderings, indistinguishable = FALSE)
#   names(fit$tcInfTransf) = names(fit$phylogenies)
#
#   fit$tcOrderings = Reduce(rbind, lapply(fit$tcInfTransf, function(w) w$drivers))
#   fit$tcOrderings = asWeightedMatrix(fit$tcOrderings)
#   fit$solutionID = solutionID
#   cat(green('OK\n'))

  # Assign what computed so far to 'x'


  # verbose = TRUE

  # Then we take all ML orderings, for every patient
  # cat(cyan(sprintf('   %20s : ', 'Expansions')), '\n')
  # x$fit$substitutions = lapply(fitPatients, transfer.orderings, x = x, verbose = verbose)
  # names(x$fit$substitutions) = fitPatients

  # print(x$fit$substitutions)

  # npatients = length(x$fit$phylogenies)
  # pb = txtProgressBar(min = 0, max = npatients, style = 3)
  # pb.status = getOption('revolver.progressBar', default = TRUE)
  #

  # x$fit$exploded = list()
  # for(patient in 1:npatients)
  # {
  #   # update progress bar
  #   if(pb.status) setTxtProgressBar(pb, patient)
  #
  #   patient = names(x$fit$phylogenies)[patient]
  #
  #   phylo = x$fit$phylogenies[[patient]]
  #   subst = x$fit$substitutions[[patient]]
  #
  #   # get groups, and the adj matrix
  #   groups = unique(phylo$dataset[phylo$dataset$is.driver, 'cluster'])
  #   adj_matrix = phylo$adj_mat
  #
  #   # visit groups in order according to a topological sort, ensures correct expansions
  #   TS = wrapTS(adj_matrix)
  #   TS = TS[TS %in% groups]
  #
  #   for(g in TS) adj_matrix = substitute.group.with.graph(g, G = adj_matrix, S = subst[[g]], verbose = verbose)
  #
  #   x$fit$exploded[[patient]] = adj_matrix
  # }
#
#
#   # Now that we have expanded, we have to update the information transfer
#   x$fit$transfer = lapply(names(x$fit$phylogenies),
#                           function(w){
#                             information.transfer_exploded_nodes(x, patient = w, transitive.orderings)
#                           } )
#   names(x$fit$transfer) = names(x$fit$phylogenies)


  class(x) = "rev_cohort_fit"

  return(x)
}




# Penalize(x, patient, E, rank = 1:3)
