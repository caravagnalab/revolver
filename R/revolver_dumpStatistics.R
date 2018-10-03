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

  pioHdr('REVOLVER: Dump statistics to Excel', "File: ", file, suffix='\t')
  pioTit('Features: drivers, occurrences & information transfer')
  ft <- revolver.featureMatrix(x)

  # helper fn for data.frame'ing
  dfWithRownames <- function(x) {
    if (is(x, "matrix")) return(cbind(Name=rownames(x), as.data.frame(x)))
    return(cbind(Name=names(x), Value=x))
  }

  sheets <- list()
  sheets[["Drivers"]] <- x$dataset[x$dataset$is.driver, -1] # wtf is Misc?
  sheets[["Occurrences"]] <- dfWithRownames(ft$occurrences)
  sheets[["Information_Transfer"]] <- ft$consensus.explosion

  if (!is.null(x$fit)) {
    pioTit('Fit: multinomial, penalty & fit')
    sheets[["Multinomial"]] <- dfWithRownames(x$fit$multinomial.penalty)
    sheets[["Penalty"]] <- dfWithRownames(x$fit$penalty)
    sheets[["Fit"]] <- dfWithRownames(Fit=x$fit$solutionID)
  }

  if (!is.null(x$cluster)) {
    pioTit('Clustering: distance, clusters & features per cluster')
    sheets[["Distance"]] <- dfWithRownames(x$cluster$distances)
    sheets[["Clusters"]] <- dfWithRownames(x$cluster$clusters)
    groups <- sort(x$cluster$clusters)
    for (g in sort(names(x$cluster$labels.colors))) {
      pts <- names(groups[groups==as.character(g)])
      cexp <- revolver.featureMatrix(x, patients=pts)$consensus.explosion
      sheets[[g]] <- dfWithRownames(cexp)
    }
  }

  if (!is.null(x$jackknife)) {
    pioTit('Jackknife: co-clustering, clusters & features per cluster')
    sheets[["Jackknife-clustering"]] <- dfWithRownames(x$jackknife$cluster)
    sheets[["Jackknife-edges-detection"]] <- dfWithRownames(x$jackknife$edges)
  }

  write_xlsx(sheets, path=file) 
  invisible(sheets)

}
