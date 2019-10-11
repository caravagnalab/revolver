#' Summary print of a cohort.
#'
#' @details Prints basic information about a cohort object,
#' status of the information available and other information.
#'
#' @param x A REVOLVER cohort.
#' @param ... Standard S3 signature.
#'
#' @return Nothing.
#' @family S3 functions
#'
#' @examples
#' # Data released in the 'evoverse.datasets'
#' data('TRACERx_NEJM_2017_REVOLVER', package = 'evoverse.datasets')
#'
#' # Cancel the fits
#' TRACERx_NEJM_2017_REVOLVER$fit = NULL
#'
#' print(TRACERx_NEJM_2017_REVOLVER)
print.rev_cohort <-
  function(x, ...)
  {
    stopifnot(inherits(x, "rev_cohort"))

    pio::pioHdr('REVOLVER - Repeated Evolution in Cancer')

    pio::pioStr("Dataset :",
                x$annotation,
                suffix = '\n',
                prefix = '\n')

    pio::pioStr(
      "Cohort  :",
      paste0(
        x$n$patients,
        ' patients, ',
        x$n$variants,
        ' variants and ',
        x$n$drivers,
        ' driver events.'
      ),
      suffix = '\n\n'
    )

    # pio::pioTit('Available computations')


    pio::pioStr('Trees per patient    :',
                ifelse(is.null(x$phylogenies), red('NO'), green('YES')),
                suffix = '\n')

    pio::pioStr('Fit via TL           :',
                ifelse(is.null(x$fit), red('NO'), green('YES')),
                suffix = '\n')

    pio::pioStr('REVOLVER clustering  :',
                ifelse(is.null(x$cluster), red('NO'), green('YES')),
                suffix = '\n')

    pio::pioStr('Jackknife statistics :',
                ifelse(is.null(x$jackknife), red('NO'), green('YES')),
                suffix = '\n')



    # pio::pioTit("Cohort summary statistics")
    # print(Stats(x))
    #
    # pio::pioTit("Drivers summary statistics")
    # print(Stats_drivers(x))


    pio::pioStr(
      "",
      "\nFor summary statistics see `?Stats_*(x)` with * = {cohort, drivers, trees, fits, clusters, ...}",
      suffix = '\n'
    )


    # if (!is.null(x$phylogenies))
    # {
    #   pio::pioTit('Summary for data and models', prefix = '\t')
    #
    #   # cat(cyan('\n\tTL Model Fit  :'), ifelse(is.null(x$fit), red('NO\n'), green('YES ')), '\n')
    #
    #   longest.name = max(nchar(names(x$phylogenies)))
    #
    #   patf = function(w, fit, clusters) {
    #     patient = names(x$phylogenies)[w]
    #
    #     s = paste(
    #       yellow(sprintf(
    #         paste('%', longest.name, 's', sep = ''), patient
    #       )),
    #       ': ',
    #       sprintf(
    #         'k = %3s | t = %3s | n = %2s | r = %2s | m = %4s | d = %2s',
    #         length(x$phylogenies[[w]]),
    #         rev_count_information_transfer_comb(x, names(x$phylogenies)[w]),
    #         x$phylogenies[[w]][[1]]$numNodes,
    #         x$phylogenies[[w]][[1]]$numRegions,
    #         nrow(x$phylogenies[[w]][[1]]$dataset),
    #         nrow(x$phylogenies[[w]][[1]]$dataset[x$phylogenies[[w]][[1]]$dataset$is.driver,])
    #       ),
    #       '\t'
    #     )
    #
    #     cat('\t', s)
    #
    #     if (fit)
    #     {
    #       stas = stats.rev_phylo(x$fit$phylogenies[[w]])
    #
    #       rank = x$fit$solutionID[w]
    #       score = format(x$fit$phylogenies[[w]]$score,
    #                      scientific = TRUE,
    #                      digits = 2)
    #       gofit = format(stas$gofit, scientific = TRUE, digits = 2)
    #
    #       cat(
    #         bgBlue('[ Fit ]', sprintf(' # %3s', yellow(rank)), ' '),
    #         '| g =',
    #         red(sprintf('%8s', gofit)),
    #         '| f =',
    #         green(sprintf('%8s', score)),
    #         '  '
    #       )
    #     }
    #
    #     if (clusters)
    #     {
    #       cluster = x$cluster$clusters[patient]
    #
    #       cat(bgMagenta('[ Cluster', sprintf('%3s', yellow(cluster)), ']'))
    #
    #     }
    #
    #     cat('\n')
    #   }
    #
    #   cat('\n')
    #   sapply(
    #     1:length(x$phylogenies),
    #     patf,
    #     fit = !is.null(x$fit),
    #     cluster = !is.null(x$cluster)
    #   )
    #
    #
    #   cat(
    #     '\n\tLegend \n\t\t k : phylogenies  \n\t\t t : combinations of information transfer \n\t\t n : groups (nodes of the tree) \n\t\t r : regions (inputs per patient) \n\t\t m : number of alterations \n\t\t d : number of driver alterations\n'
    #   )
    #
    #   if (!is.null(x$fit))
    #     cat(
    #       '\t\t # : number of the solution selection (out of k)  \n\t\t g : goodness-of-fit \n\t\t f : score of the model\n'
    #     )
    #
    #   if (!is.null(x$cluster))
    #   {
    #     pio::pioTit('Summary for clusters', prefix = '\t')
    #
    #     ids = names(x$cluster$labels.colors)
    #
    #     medians = rep(NA, length(ids))
    #
    #     if (!is.null(x$jackknife)) {
    #       medians = x$jackknife$cluster.medians
    #     }
    #
    #     size = as.vector(table(x$cluster$clusters))
    #
    #     df = data.frame(Cluster = ids,
    #                     n = size,
    #                     Jackknife = medians)
    #     apply(df, 1, function(w) {
    #       cat('\n\t',
    #           paste(
    #             yellow(sprintf('%4s', w['Cluster'])),
    #             ': ',
    #             sprintf('n = %4s | jackknife = %4s',
    #                     w['n'],
    #                     w['Jackknife'])
    #           ))
    #
    #     })
    #
    #     cat('\n\n\tLegend \n\t\t         n : number of patients in the cluster')
    #
    #     if (!is.null(x$jackknife))
    #       cat('\n\t\t jackknife : median co-clustering probability estiamted with jackknife\n')
    #
    #     # assignments = x$cluster$clusters
    #     # edges = x$jackknife$edges
    #     # pio::pioTit('Summary for edges', prefix = '\t')
    #     # colnames(edges)[3] = 'Jackknife'
    #     # pio::pioDisp(edges)
    #   }
    # }

    cat('\n')
    revolver_check_cohort(x)

  }

#' Summary print of a cohort with fits.
#'
#' @details Prints basic information about a cohort object,
#' status of the information available and other information.
#'
#' @param x A REVOLVER cohort with fits.
#' @param ... Standard S3 signature.
#'
#' @return Nothing
#' @family S3 functions
#'
#' @examples
#' # Data released in the 'evoverse.datasets'
#' data('TRACERx_NEJM_2017_REVOLVER', package = 'evoverse.datasets')
#'
#' print(TRACERx_NEJM_2017_REVOLVER)
print.rev_cohort_fit = function(x, ...)
{
  class(x) = 'rev_cohort'
  print(x)
}


