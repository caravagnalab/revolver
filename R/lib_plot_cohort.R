## Plots for a cohort are of different kinds.


#' @title Plot a patient's data
#'
#' @param x An object of class \code{"rev_cohort"}
#' @param patient The patient to plot
#' @param cex Cex for graphics (pheatmap cells)
#' @param file Output file, if not NA.
#'
#' @return nothing
#' @export
#' @import crayon
#'
#' @examples
#' data(Breast.fit)
#' revolver_plt_patient_data(Breast.fit, 'PD14767')
revolver_plt_patient_data = function(x,
                                      patient,
                                      cex = 1,
                                      file = NA)
{
  pio::pioHdr(paste('REVOLVER Plot: data for', patient), toPrint = NULL)

  pat = x$phylogenies[[patient]][[1]]
  ccf = pat$CCF[, pat$samples_ID, drop = F]

  ccf.numbers = ccf
  ccf.numbers = round(ccf.numbers, 3)
  ccf.numbers[ccf.numbers == 0] = ''

  annot = pat$CCF[, c('is.driver', 'is.clonal')]
  annot$is.clonal = as.character(annot$is.clonal)
  annot$is.driver = as.character(annot$is.driver)

  yes_no_col = c('forestgreen', 'gainsboro')
  names(yes_no_col) = c('TRUE', 'FALSE')
  ann_colors = list(is.driver = yes_no_col, is.clonal = yes_no_col)

  br = c(0, 1e-3, seq(0.1, 1, 0.1))

  pheatmap::pheatmap(
    ccf,
    main = paste('Data for', patient),
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    border = NA,
    breaks = br,
    legend = F,
    color = c('white', scols(1:(length(
      br
    ) - 1), 'Blues')),
    display_numbers = ccf.numbers,
    cellwidth = 34 * cex,
    cellheight = 12 * cex,
    number_color = 'orange',
    annotation_row = annot,
    annotation_colors = ann_colors,
    file = file
  )
}



#' Plot the scores for the trees of each patient.
#'
#' @param x An object of class \code{"rev_cohort"}
#' @param patient The patient to plot
#' @param cex Cex for graphics (pheatmap cells)
#' @param file Output file, if not NA.
#'
#' @return nothing
#' @export
#' 
#' @importFrom graphics barplot
#'
#' @examples
#' data(Breast.fit)
#' revolver_plt_patient_trees_scores(Breast.fit, 'PD14767')
revolver_plt_patient_trees_scores = function(x,
                                             patient,
                                             cex = 1,
                                             file = NA)
{
  obj_has_trees(x)
  pio::pioHdr(paste('REVOLVER Plot: tree scores for', patient), toPrint = NULL)

  scores = lapply(x$phylogenies[[patient]], function(e)
    e$score)
  comb = rev_count_information_transfer_comb(x, patient)
  colors = scols(1:comb, 'Dark2')

  groups = lapply(x$phylogenies[[patient]],
                  function(w)
                    paste(sort(DataFrameToEdges(
                      w$transfer$drivers
                    )), collapse = ':'))
  groups = unlist(groups)
  names(colors) = unique(groups)

  lot = mylayout.on(file, 1, c(4, 2), cex)

  barplot(
    unlist(scores),
    col = colors[groups],
    border = NA,
    main = paste("Scored trees for", patient),
    sub = paste("Number of combinations of transfer", comb),
    horiz = T,
    ylab = 'Trees',
    xlab = 'Score'
  )

  # legend('topright', legend = comb, bty = 'n')

  mylayout.off(lot)
  invisible(NULL)
}

#' @title Plot a the top-scoring trees for a patient.
#'
#' @details Iterative plotting for the set of input trees of a patient (mutation
#' trees or phylogenetic trees). This is usually used when ispecting the
#' information transfer for the input patients.
#'
#' @param x An object of class \code{"rev_cohort"}
#' @param patient The patient to plot.
#' @param max.phylogenies How many trees should be computed.
#' @param file Output file, or NA.
#' @param cex Scale cex for graphics.
#'
#' @return nothing
#' @export
#'
#' @examples
#' data(Breast.fit)
#' revolver_plt_patient_trees(Breast.fit, 'PD14767')
revolver_plt_patient_trees = function(x,
                                      patient,
                                      max.phylogenies = 12,
                                      file = NA,
                                      cex = 1)
{
  obj_has_trees(x)

  ##### Plot the distribution
  distribution = x$phylogenies[[patient]]
  ndistribution = length(distribution)

  if (ndistribution > max.phylogenies)
  {
    distribution = distribution[1:max.phylogenies]
  }
  else
  {
    max.phylogenies = ndistribution
  }

  pio::pioHdr(paste('REVOLVER Plot: top', max.phylogenies, 'ranked trees for', patient, '-- out of', ndistribution), toPrint = NULL)
  cat("Progress: ")

  lot = mylayout.on(file, length(distribution), c(4, 4), cex)

  toMerge = sapply(1:length(distribution),
                   function(x) {
                      plot(
                       distribution[[x]],
                       # plot.stat = plot.stat,
                       cex = cex
                       # table.cex = cex
                       # edge.label = round(MI.table, 2),
                     )

                     cat(paste('.', x, sep  = ''))
                   })

  mylayout.off(lot)

  cat('\n')

  invisible(NULL)
}
