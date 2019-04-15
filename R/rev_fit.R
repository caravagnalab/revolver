
#' Print a REVOLVER cohort object with fits.
#'
#' @param x A \code{"rev_cohort_fit"} object
#' @param ... unused arguments from generic
#'
#' @return none
#' @export print.rev_cohort_fit
#'
#' @examples
#' data(Breast.fit)
#' Breast.fit
print.rev_cohort_fit = function(x, ...)
{
  class(x) = 'rev_cohort'
  print(x)
}

#' @title Plot fits from a REVOLVER cohort.
#'
#' @details
#' Iterative plotting functions that scans a cohort and runs \code{\link{revolver_report_fit_patient}}
#' on a set of patients (default all).
#'
#' @param x An object of class \code{"rev_cohort_fit"}
#' @param patients The patients to plot, default is all the one available.
#' @param merge.files Merge the output files? (FALSE)
#' @param cex Scale cex for graphics
#'
#' @return nothing
#' @export plot.rev_cohort_fit
#' @import crayon
#'
#' @examples
#' data(Breast.fit)
#' plot.rev_cohort_fit(Breast.fit, patients = Breast.fit$patients[1:4])
plot.rev_cohort_fit = function(x,
                           patients = x$patients,
                           merge.files = FALSE,
                           cex = 1)
{
  obj_has_trees(x)

  pio::pioHdr('REVOLVER Plot: fits (tree, trajectories, information transfer)',
              c(
                `Patients` = paste(patients, collapse = ', '),
                `Merge output PDFs` = merge.files
              ),
              prefix = '\t -')

  # if(is.na(file)) stop('A file is required for this plot!')

  for (patient in patients)
  {
    pio::pioTit(paste("Processing", patient))

    revolver_report_fit_patient(x, patient = patient, cex = cex,
                                file = paste0('REVOLVER-report-fit-patient-', patient, '.pdf'))
  }

  if(merge.files)
    jamPDF(
      in.files = paste0('REVOLVER-report-fit-patient-', patients, '.pdf'),
      out.file = paste0('REVOLVER-report-fit-all-cohort-merged.pdf'),
      delete.original = FALSE,
      layout = '1x1'
      )


  invisible(NULL)
}



#' @title REVOLVER fit function.
#'
#' @details  The main fitting function for REVOLVER, implements the 2-steps EM-algorithm described
#' in the main paper.
#'
#' @param x A \code{"rev_cohort"} object with available trees.
#' @param initial.solution Either a scalar to fix one initial condition (rank id), or \code{NA} to sample it randomly.
#' @param max.iterations Maximum number of EM steps before forcing stop.
#' @param restarts Number of initial conditions sampled to compute optimal fit
#' @param parallel Use a parallel engine with \code{TRUE}, hides ouptut. Exploits \code{doParallel}.
#' @param cores.ratio Percentage of the number of cores to use, after the overall number of
#' available cores has been detected by \code{doParallel}.
#' @param transitive.orderings If transitive orderings should be computed before the second
#' fitting step (transfering of orderings). The default is \code{FALSE}, which we suggest.
#' @param verbose Verbose output
#'
#' @return A \code{"rev_cohort_fit"} object
#' @export
#' @import parallel
#' @import crayon
#' @import doParallel
#' @import foreach
#'
#'
#' @examples
#' data(CRC.cohort)
#' fit = revolver_fit(CRC.cohort)
# revolver_fit = function(x,
#                         initial.solution = 1,
#                         max.iterations = 10,
#                         restarts = 0,
#                         parallel = FALSE,
#                         cores.ratio = .8,
#                         transitive.orderings = FALSE,
#                         verbose = FALSE)
# {
#   # Fancy prints
#   pline = function() {cat('\n\n\t'); for(i in 1:20)cat((cyan('~/-\\~')))}
#
#   # Program name
#   cat(bgCyan('REVOLVER'), cyan(x$REVOLVER_VERSION))
#   cat(blue(paste(' -- ', x$annotation, sep ='')), '\n')
#
#   revolver_check_cohort(x, auto.fix = FALSE)
#
#   # Title
#   cat(yellow('\n\t== Fitting via Transfer Learning == \n'))
#   if(verbose) cat(yellow('\n\t== VERBOSE is ON == \n'))
#
#   numPatients = length(x$phylogenies)
#   fitPatients = names(x$phylogenies)
#
#   # Input patients
#   cat(cyan('\nFitting '), 'n =', length(x$patients), cyan('patients \n\n   '))
#   sapply(1:length(fitPatients), function(p) {
#     cat(sprintf('%20s', fitPatients[p]))
#     if(p %% 5 == 0) cat('\n   ')
#   })
#
#   # Initial condition
#   cat(cyan('\nInitial solution :'), ifelse(!is.na(initial.solution), 'Fixed', 'Random\n'))
#
#   # Restarts
#   if(restarts > 0 & !is.na(initial.solution)) {
#     cat(red(' -- You have set an initial condition, restarts and/ or paralllel will be disregarded.'))
#     restarts = 0
#     paralllel = FALSE
#   }
#
#   if(restarts > 0)
#   {
#     cat(cyan('\nIID initial conditions :'), restarts, cyan('\t parallel'),
#         ifelse(parallel, green('TRUE'), red('FALSE')), '\n')
#   }
#
#   results = NULL
#
#   # Parallel fit
#   if(parallel) {
#
#     pclusters = setup_parallel(cores.ratio)
#
#     r = foreach(num = 0:restarts, .packages = c("crayon", "igraph", 'matrixStats', 'matrixcalc'), .export = ls(globalenv())) %dopar%
#     {
#       tl_revolver_fit(x, initial.solution, max.iterations, transitive.orderings, verbose)
#     }
#
#     results = r
#
#     stop_parallel(pclusters)
#   }
#   else {
#     # NON Parallel fit
#
#     for(i in 0:restarts) {
#       pline()
#       model_fit = tl_revolver_fit(x, initial.solution, max.iterations, transitive.orderings, verbose)
#       results = append(results, list(model_fit))
#     }
#   }
#
#   # Get the best solution
#   if(restarts > 0)
#   {
#     pline()
#     cat(cyan('\n\nMedian goodness-of-fit (from fit penalty): '))
#     scores = sapply(results, function(w) median(w$fit$penalty))
#     cat(paste(scores, collapse = ', '))
#     cat(cyan('\t Best is'), bgGreen('',scores[which.max(scores)], ''), '\n')
#     return(results[[which.max(scores)]])
#   }
#   else return(results[[1]])
# }



ML.orderings = function(group, orderings)
{
  drivers.list = group$variantID
  M = orderings[drivers.list, drivers.list, drop = FALSE]
  M.all = M

  ml.parents = function(x, n){
    par = pi(x, n)
    par = which(x[, n] < max(x[, n]))
    if(length(par) > 0) x[par, n] = 0
    return(x)
  }

  if(length(drivers.list) > 1 && !all(M == 0))
    for(d in drivers.list) M = ml.parents(M, d)

  return(list(
    ML = M,
    all = M.all
  ))
}

# Transfer the ordering from other patients to this one
transfer.orderings = function(x, patient, verbose = FALSE)
{
  patient.id = patient
  patient.phylo = x$fit$phylogenies[[patient]]
  drivers.list = patient.phylo$dataset[patient.phylo$dataset$is.driver, ]
  drivers.list = drivers.list[order(drivers.list$cluster), ]

  # print(drivers.list)
  # print(split(drivers.list, f = drivers.list$cluster))

  # every group is expanded independently
  ret = lapply(
    split(drivers.list, f = drivers.list$cluster),
    function(group, orderings)
    {
      # The current group is a cluster
      cluster.id = as.character(unique(group$cluster))
      ord = ML.orderings(group, orderings)

      M = ord$ML
      M.all = ord$all

      # remove this so far...
      # TS = topo_sort(M, mode = 'out')

      if(verbose) {
        GM = graph_from_adjacency_matrix(M, weighted = TRUE)
        GM.all = graph_from_adjacency_matrix(M.all, weighted = TRUE)

        cat(cyan('\nML parents:\n'))
        print(M)

        # cat(cyan('\nTopological sort:'), paste(TS), ':',  paste(V(M)$name[TS], collapse = ' --> '), '\n')

        pdf('tmp.pdf', width = 5, height = 5)

        plot(GM,
             layout = layout_in_circle(GM),
             vertex.size = 20,
             vertex.color = 'steelblue',
             vertex.frame.color = 'white',
             edge.label = igraph::E(GM)$weight)
        title(paste('Group ', cluster.id,  ' - Mult. counts'))

        plot(GM.all,
             layout = layout_in_circle(GM.all),
             vertex.size = 20,
             vertex.color = 'steelblue',
             vertex.frame.color = 'white',
             edge.label = igraph::E(GM.all)$weight)
        title('Max. Lik. parents')
        #
        # TSe = as.matrix(as_adjacency_matrix(M))
        # TSe[T] = 0
        #
        # if(sum(as_adjacency_matrix(M)) > 0) for(i in 2:length(TS)) TSe[TS[i-1], TS[i]] = 1
        #
        # plot(graph_from_adjacency_matrix(TSe),
        #      vertex.color = 'steelblue',
        #      vertex.frame.color = 'white'
        # )
        # title('Topological sort ')

        dev.off()

        jamPDF('tmp.pdf', paste(cluster.id, '.pdf', sep = ''), layout = '2x1')
      }

      return(M)
    },
    orderings = x$fit$tcOrdering)

  # print(ret)

  names(ret) = unique(drivers.list$cluster)

  if(verbose)
    jamPDF(
      paste(unique(drivers.list$cluster), '.pdf', sep = ''),
      paste(patient, '.pdf', sep = ''), layout = '1x1')

  return(ret)
}

optimize = function(E, solutions, current.one, alpha = 1, greedy = FALSE, verbose = FALSE)
{
  add.to.Matrix = function(M, df){
    for(i in 1:nrow(df))
      M[df[i, 'from'], df[i, 'to']] = M[df[i, 'from'], df[i, 'to']] + 1
    M
  }

  del.from.Matrix = function(M, df){
    for(i in 1:nrow(df))
      M[df[i, 'from'], df[i, 'to']] = M[df[i, 'from'], df[i, 'to']] - 1
    M
  }

  penalize = function(transfer, M) {
    cbind(
      transfer,
      penalty = apply(transfer,
                      1,
                      function(e){
                        ifelse(
                          M[e[1], e[2]] == 0,
                          1,
                          M[e[1], e[2]]
                        )/ sum(M[, e[2]]) }
      ))
    }
  ######


  current.transfer = penalize(solutions[[current.one]]$transfer$drivers, E)
  current.score = solutions[[current.one]]$score * prod(current.transfer$penalty)^alpha
  E.cached = del.from.Matrix(E, current.transfer)
  #
  # print(penalize(solutions[[current.one]]$transfer$drivers, E))
  # print(solutions[[current.one]]$transfer$drivers)
  # print(E.cached[,'MHK21'])
  #

  if(verbose)
  {
    cat(cyan('Optimizer for model :'), yellow(solutions[[current.one]]$patient_ID), cyan('\t Greedy'), ifelse(greedy, green(TRUE), red(FALSE)), '\n',
        cyan('Solution ='), current.one, ' out of', length(solutions), '\n',
        cyan('       f ='), blue(current.score), '\n',
        cyan('   Score ='), solutions[[current.one]]$score, '\n',
        cyan(' Penalty ='), prod(current.transfer$penalty), '\n')
    print(current.transfer)
  }

  best = c(current.score, current.one, prod(current.transfer$penalty)^alpha)
  names(best) = c('score', 'idx', 'penalty')

  # ML estimate of a single patient, given the E-matrix E
  for(i in 1:length(solutions))
  {
    if(i == current.one) next;

    # this is what we test
    solution = solutions[[i]]

    # first, we compute the differential counts for E if we had this solution, instead of the current
    E.new = add.to.Matrix(E.cached, solution$transfer$drivers)

    penalty.new = penalize(solution$transfer$drivers, E.new)
    f.new = prod(penalty.new$penalty)^alpha * solution$score

    # print(penalty.new)

    if(verbose) {
      cat(inverse('\nTesting solution #'),  yellow(sprintf('%4s', i)),
          cyan(sprintf('%10s', 'Score')), sprintf('%10s', round(solution$score, 6)),
          cyan(sprintf('%10s', 'Penalty')), sprintf('%10s', round(prod(penalty.new$penalty)^alpha, 6)),
          cyan(sprintf('%10s', 'f')), ifelse(f.new > best['score'], green(f.new), red(f.new))
          )

      # print(penalty)

    }

    if(f.new > best['score'])
    {
      best['score'] = f.new
      best['idx'] = i
      best['penalty'] = prod(penalty.new$penalty)^alpha

      if(greedy) break;
    }

  }

  return(best)
}

rev_table.orderings = function(x, intelli.cutoff = 3)
{
  if(is.null(x$fit)) stop('No fit in this object?')

  where.is = function(from, to) {
   idx = sapply(edges.conv, function(w) paste(from, to, sep = '~') %in% w)
   samples = sort(names(edges.conv)[idx])

   out.samples = NULL
   counts = 0
   for(s in samples) {
     drv = x$fit$phylogenies[[s]]$transfer$drivers

      if(!(paste(from, to, sep = '~') %in% DataFrameToEdges(drv))) {
        counts = counts + 1
      } else  out.samples = c(out.samples, s)

   }

   list(out.samples, counts)
  }

  edges.conv = lapply(x$fit$transfer, DataFrameToEdges)
  edges.preTr = lapply(x$fit$phylogenies, function(w) DataFrameToEdges(w$transfer$drivers))

  edges = unlist(edges.conv)
  edges = table(edges)

  tab = Reduce(rbind, x$fit$transfer)
  tab = unique(tab)

  keys = DataFrameToEdges(tab)
  tab$Count = edges[keys]
  tab = tab[order(tab$Count, decreasing = TRUE), , drop = FALSE]
  rownames(tab) = NULL

  intelli = tab[tab$Count >= intelli.cutoff, , drop = FALSE]
  intelli = split(intelli, f = intelli$to)
  intelli = lapply(intelli, function(w){
    n = apply(w, 1, function(z)
      {
      wh = where.is(z['from'], z['to'])
      data.frame(
        Count.expanded = wh[[2]],
        patientID = paste(wh[[1]], collapse = ', ')
        )
      })
    n = Reduce(rbind, n)

    w = cbind(w, n)
    return(w)
  })

  #
  # intelli.filter = lapply(intelli, function(w){
  #   if(nrow(w) == 1) return(TRUE)
  #   max = w[which.max(w$Count)[1], ]
  #   others = w[setdiff(1:nrow(w), which.max(w$Count)), ]
  #   return(max$Count > sum(others$Count))
  # })

  return(list(tab = tab, intelli = intelli))
}

