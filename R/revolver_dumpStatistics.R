#' Dump statistics for fit to an Excel file.
#'
#' @param x REVOLVER cohort fit
#' @param file output file
#'
#' @return invisibly, the list of data.frames written to the Excel file 
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

  sheets <- list()
  pioHdr('REVOLVER: Dump statistics to Excel', "File: ", file, suffix='\t')
  pioTit('Features: drivers, occurrences & information transfer')
  ft <- revolver.featureMatrix(x)

  sheets[["Drivers"]] <- x$dataset[x$dataset$is.driver,]
  sheets[["Occurrences"]] <- as.data.frame(ft$occurrences)
  sheets[["Information_Transfer"]] <- ft$consensus.explosion

  if (!is.null(x$fit)) {
    pioTit('Fit: multinomial, penalty & fit')
    sheets[["Multinomial"]] <- as.data.frame(x$fit$multinomial.penalty)
    sheets[["Penalty"]] <- data.frame(Penalty=x$fit$penalty)
    sheets[["Fit"]] <- data.frame(Fit=x$fit$solutionID)
  }

  if (!is.null(x$cluster)) {
    pioTit('Clustering: distance, clusters & features per cluster')
    sheets[["Distance"]] <- data.frame(Distance=x$cluster$distances)
    sheets[["Clusters"]] <- data.frame(Clusters=x$cluster$clusters)
    groups <- sort(x$cluster$clusters)
    for (g in sort(names(x$cluster$labels.colors))) {
      pts <- names(groups[groups==as.character(g)])
      cexp <- revolver.featureMatrix(x, patients=pts)$consensus.explosion
      sheets[[g]] <- as.data.frame(cexp)
    }
  }

  if (!is.null(x$jackknife)) {
    pioTit('Jackknife: co-clustering, clusters & features per cluster')
    sheets[["Jackknife-clustering"]] <- data.frame(Cluster=x$jackknife$cluster)
    sheets[["Jackknife-edges-detection"]] <- data.frame(Edges=x$jackknife$edges)
  }

  write_xlsx(sheets, path=file) 
  invisible(sheets)

}
