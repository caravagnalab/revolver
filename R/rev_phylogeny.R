



#' @title Summary of a \code{"rev_phylo"} object.
#' @details  Reports some summary for a \code{"rev_phylo"} object, this is a bit more than just using print.
#'
#' @param x An object of class \code{"rev_phylo"}.
#' @param ... Extra parameters
#' @param print.stat Print or not statistics to screen.
#'
#' @return none
#' @export summary.rev_phylo
#' @import crayon
#'
#' @examples
#' data(CRC.cohort)
#' summary(CRC.cohort$phylogenies[['adenoma_3']][[1]])
summary.rev_phylo <- function(x, ..., print.stat = TRUE) {
  print.rev_phylo(x, ...)

  cat(cyan('\nCCF clusters:\n'))
  print(x$CCF)

  cat(cyan('\nDrivers:\n'))
  print(x$dataset[x$dataset$is.driver,])

  cat(cyan('\nInformation Transfer:\n'))
  print(x$transfer)

  stat = stats.rev_phylo(x)

  if (print.stat)  {
    cat(yellow('\n[print.stat = T] Empirical probabilities:\n'))
    print.noquote(stat$probs)

    cat(yellow('\n[print.stat = T] Suppes\' conditions (binary data):\n'))
    print.noquote(stat$Suppes)

    # cat(yellow('\n[print.stat = T] Crossing rule (CCF):\n'))
    # print.noquote(stat$CCF.crossing.rule)

    cat(yellow('\n[print.stat = T] Pigeonhole principle (CCF):\n'))
    print.noquote(stat$CCF.pigeonhole)
  }

  cat(cyan('\nViolations:\n'))
  cat(
    'Suppes\' Temporal Priority       : ',
    green(nrow(stat$Suppes) - stat$violations['tp']),
    red(stat$violations['tp']),
    '\n'
  )
  # cat('\t\tCCF Crossing Rule       : ', green(nrow(stat$CCF.crossing.rule) - stat$violations['cr']), red(stat$violations['cr']), '\n')
  cat(
    'Suppes\' Probability Raising     : ',
    green(nrow(stat$Suppes) - stat$violations['pr']),
    red(stat$violations['pr']),
    '\n'
  )
  cat(
    'CCF Pigenhole Principle         : ',
    green(
      nrow(stat$CCF.pigeonhole) * ncol(stat$CCF.pigeonhole) - stat$violations['pp']
    ),
    red(stat$violations['pp']),
    '\n'
  )

  cat(cyan('\nGoodness-of-fit: '), stat$gofit)
}





#' @title Compute statistics for a \code{"rev_phylo"} tree.
#'
#' @details
#' Compute some statistics for a \code{"rev_phylo"} tree. Thee include violations of the
#' pigeonhole principle, Suppes' conditions etc. Stats are compute no matter what
#' data is used to create the trees, but one should only consider those relevant
#' to the actual type of data.
#'
#' @param x An object of class \code{"rev_phylo"}.
#'
#' @return A list with all computed statistics.
#' @export
#'
#' @examples
#' data(CRC.cohort)
#' stats.rev_phylo(CRC.cohort$phylogenies[['adenoma_3']][[1]])
stats.rev_phylo = function(x)
{
  probs = function(b, i, j) {
    b = rbind(b, GL = 0)
    mi = sum(b[, i]) / nrow(b)
    mj = sum(b[, j]) / nrow(b)
    joi = (b[, i] %*% b[, j]) / nrow(b)
    return(c(mi = mi, mj = mj, joi = joi))
  }

  M = MatrixToDataFrame(x$adj_mat)

  if (x$germline)
    M = M[M$from != 'GL', , drop = FALSE]

  if (nrow(M) == 0) {
    return(
      list(
        probs = NULL,
        violations = NULL,
        gofit = 1,
        Suppes = NULL,
        CCF.crossing.rule = NULL,
        CCF.pigeonhole = NULL
      )
    )

  }

  p = data.frame(matrix(nrow = 0, ncol = 3), stringsAsFactors = FALSE)
  colnames(p) = c('p(i)', 'p(j)', 'p(i,j)')

  p = apply(M, 1,
            function(e) {
              probs(x$binary.data, e[1], e[2])
            })
  p = t(p)
  rownames(p) = apply(M, 1, function(e)
    paste(e[1], '->', e[2]))
  colnames(p) = c('p(i)', 'p(j)', 'p(i,j)')

  TP = apply(p,
             1,
             function(e) {
               if (e[1] > e[2])
                 return('>')
               if (e[1] == e[2])
                 return('=')
               return('<')
             })

  PR = apply(p,
             1,
             function(e) {
               if (e[1] * e[2] < e[3])
                 return('>')
               if (e[1] * e[2] == e[3])
                 return('=')
               return('<')
             })

  df.Suppes = cbind(TP, PR)
  colnames(df.Suppes) = c('Temporal Priority', 'Probability Raising')

  # CCF data for this patient
  CCF = x$CCF[, x$samples_ID, drop = FALSE]

  # Direction penalty -- it was embedded in PP so it's now commented
  # direction.penalty = apply(
  #   M,
  #   1,
  #   function(e){
  #   CCF[e[1], , drop = FALSE] >= CCF[e[2], , drop = FALSE]
  # })
  #
  # if(ncol(CCF) == 1) direction.penalty = matrix(direction.penalty, ncol = 1)
  # else   direction.penalty = t(direction.penalty)
  # rownames(direction.penalty) = rownames(p)
  # colnames(direction.penalty) = x$samples_ID

  ## Branching penalty is the Pigeonhole Principle
  Mmatrix = DataFrameToMatrix(M)

  branching.penalty = sapply(x$nodes_ID,
                             function(n) {
                               # print(n)
                               cl = children(Mmatrix, n)
                               if (length(cl) == 0)
                                 return(NULL)
                               CCF[n, , drop = FALSE] >= colSums(CCF[cl, , drop = FALSE])
                             })
  branching.penalty = branching.penalty[!unlist(lapply(branching.penalty, is.null))]
  branchp = Reduce(rbind, branching.penalty)
  rownames(branchp) = names(branching.penalty)

  # Scores summary -- we removed the direction penalty
  tp.viol = sum(as.integer(df.Suppes[, 1] == '<'))
  pr.viol = sum(as.integer(df.Suppes[, 2] == '<'))
  # cr.viol = sum(as.integer(!direction.penalty))
  ph.viol = sum(as.integer(!branchp))

  # violations = c(tp = tp.viol, pr = pr.viol, cr = cr.viol, pp = ph.viol)
  violations = c(tp = tp.viol, pr = pr.viol, pp = ph.viol)

  # Some alternative Goodness-Of-Fit measures. In the end, we take the
  # exact definition of the penalty for phylogenetic trees.
  # gofit = 2 * (x$numNodes - 1) * x$numRegions + x$numNodes  * x$numRegions
  # gofit = 1 - sum(violations)/gofit
  # gofit = nrow(branchp) * x$numRegions
  # gofit = 1 - ph.viol/gofit
  capture.output({
    gofit = node.penalty.for.branching(list(M), CCF)
  })


  return(
    list(
      probs = p,
      violations = violations,
      gofit = gofit,
      Suppes = df.Suppes,
      # CCF.crossing.rule = direction.penalty, # Removed
      CCF.pigeonhole = branchp
    )
  )

}



#'
#' #' S3 method calling \code{\link{revolver_plt_tree}}.
#' #'
#' #' @param x An object of class \code{"rev_phylo"}.
#' #' @param file Output file, or \code{NA}.
#' #' @param palette RColorBrewer palette to colour clusters.
#' #' @param cex Cex for the graph.
#' #' @param alpha Transparency.
#' #' @param verbose Output.
#' #' @param ... Extra parameters
#' #'
#' #' @return Nothing
#' #' @export plot.rev_phylo
#' #' @import crayon
#' #'
#' #' @examples
#' #' data(CRC.cohort)
#' #' plot(CRC.cohort$phylogenies[['adenoma_3']][[1]])
#' plot.rev_phylo = function(x,
#'                           file = NA,
#'                           palette = 'Set1',
#'                           cex = 1,
#'                           alpha = 0.7)
#' {
#'   revolver_plt_tree(
#'     x,
#'     file = file,
#'     # edge.width = edge.width,
#'     # edge.label = edge.label,
#'     palette = palette,
#'     cex = cex,
#'     alpha = alpha
#'   )
#'   invisible(NULL)
#' }

# # a-la-dictionary
# driver = function(x, c) {
#   x$dataset[which(x$dataset$cluster == c &
#                     x$dataset$is.driver), 'variantID']
# }
# driver.indexOf = function(x, d) {
#   x$dataset[which(x$dataset$variantID == d &
#                     x$dataset$is.driver), 'cluster']
# }

# # Compute the frontier of "var" in a model. The frontier is the set of
# # mutations in Gamma that are directly reachable via
# # the transitive closure of the --> relation. So, these are the events
# # selected by "var"'s evolutionary trajectories.
# information.transfer = function(x,
#                                 transitive.closure = FALSE,
#                                 indistinguishable = FALSE)
# {
#   aux = function(r)
#   {
#     r.d = driver(x, r)
#
#     if (length(r.d) > 0)
#       return(r) # stop recursion
#
#     c = children(model, r)
#
#     if (is.null(c))
#       return(NULL) # leaf
#
#     # recursion, reduction etc.
#     return(Reduce(union,
#                   sapply(c, aux)))
#   }
#
#   model = x$adj_mat
#   nodes.drivers = x$dataset %>%
#     filter(is.driver) %>%
#     pull(cluster)
#
#   # Root
#   df = expand.grid(
#     from = x$root,
#     to = aux(x$root),
#     stringsAsFactors = FALSE
#   )
#
#   # and all other stuff
#   for (n in nodes.drivers)
#   {
#     df.n = NULL
#
#     for (c in children(model, n))
#       df.n = rbind(df.n,
#                    expand.grid(
#                      from = n,
#                      to = aux(c),
#                      stringsAsFactors = FALSE
#                    ))
#
#     df = rbind(df, df.n)
#   }
#
#   # Then, for all clones, we expand the drivers that they have
#   expanded = apply(df,
#                    1,
#                    function(w) {
#                      e.from = driver(x, w['from'])
#                      if (length(e.from) == 0)
#                        e.from = w['from']
#
#                      e.to = driver(x, w['to'])
#
#                      expand.grid(from = e.from,
#                                  to = e.to,
#                                  stringsAsFactors = FALSE)
#                    })
#   expanded = Reduce(rbind, expanded)
#
#   # Then we see if we want to expand also the indistinguishible
#   if (indistinguishable) {
#     for (n in nodes.drivers) {
#       expanded = rbind(expanded,
#                        expand.grid(
#                          from = driver(x, n),
#                          to = driver(x, n),
#                          stringsAsFactors = FALSE
#                        ))
#     }
#   }
#
#   # If we need transitive closures, we compute them here
#   if (transitive.closure)
#   {
#     df = DataFrameToMatrix(df)
#     df = nem::transitive.closure(df, mat = T, loops = FALSE)
#
#     expanded = DataFrameToMatrix(expanded)
#     expanded = nem::transitive.closure(expanded, mat = T, loops = FALSE)
#
#     df = MatrixToDataFrame(df)
#     expanded = MatrixToDataFrame(expanded)
#   }
#
#   return(list(clones = df, drivers = expanded))
# }


