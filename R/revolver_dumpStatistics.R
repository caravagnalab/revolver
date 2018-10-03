#' Dump statistics for fit to an Excel file.
#'
#' @param x REVOLVER cohort fit
#' @param file output file
#'
#' @return none
#' 
#' @importFrom writexl write_xlsx
#' @importFrom pio pioHdr pioTit
#'
#' @examples
#' data(CRC.cohort)
#' fit <- revolver_fit(CRC.cohort)
#' revolver_dumpStatistics(fit)
#'
#' @export
revolver_dumpStatistics <- function(x, file='REVOLVER-Statistics.xlsx') {

  pioHdr('REVOLVER: Dump statistics to Excel', "File: ", file, suffix='\t')
  pioTit('Features: drivers, occurrences & information transfer')
  ft <- revolver.featureMatrix(x)

  sheets <- list(
   "Drivers"=x$dataset[x$dataset$is.driver,],
   "Occurrences"=ft$occurrences,
   "Information_Transfer"=ft$consensus.explosion,
  )

  if (!is.null(x$fit)) {
    pioTit('Fit: multinomial, penalty & fit')
    sheets[["Multinomial"]] <- x$fit$multinomial.penalty
    sheets[["Penalty"]] <- x$fit$penalty
    sheets[["Fit"]] <- x$fit$solutionID
  }

  if (!is.null(x$cluster)) {
    pioTit('Clustering: distance, clusters & features per cluster')
    sheets[["Distance"]] <- x$cluster$distances
    sheets[["Clusters"]] <- x$cluster$clusters
    groups <- sort(x$cluster$clusters)
    for (g in sort(names(x$cluster$labels.colors))) {
      pts <- names(groups[groups==as.character(g)])
      sheets[[g]] <- revolver.featureMatrix(x, patients=pts)$consensus.explosion
    }
  }

  if (!is.null(x$jackknife)) {
    pioTit('Jackknife: co-clustering, clusters & features per cluster')
    sheets[["Jackknife-clustering"]] <- x$jackknife$cluster
    sheets[["Jackknife-edges-detection"]] <- x$jackknife$edges
  }

  write_xlsx(sheets, path=file) 

}

