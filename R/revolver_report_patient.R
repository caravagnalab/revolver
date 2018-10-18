#' @title Generate PDF report for a patient's data and trees
#'
#' @param x An object of class \code{"rev_cohort"}
#' @param patient The patient to plot
#' @param cex Cex for graphics
#' @param max.phylogenies Max number of trees to plot.
#' @param file Output file, cannot be NA.
#'
#' @return nothing
#' @export
#' @import crayon
#'
#' @examples
#' data(Breast.fit)
#' revolver_report_patient(Breast.fit, 'PD9770')
revolver_report_patient = function(x,
                                   patient,
                                   cex = 1,
                                   max.phylogenies = 12,
                                   file = paste0('REVOLVER-report-patient-data-models-', patient, '.pdf'))
{
  obj_has_trees(x)
  pio::pioHdr(paste('REVOLVER Plot: report on data and trees for', patient), toPrint = NULL)

  if(is.na(file)) stop("Output file cannot be NA")

  f = paste(patient,
            c('data.pdf', 'trees_scores.pdf', 'trees.pdf'), sep = '-')

  # dumping everything to a PDF all the time does NOT make debugging easier... 
  revolver_plt_patient_data(x, patient, cex = cex, file = f[1])
  revolver_plt_patient_trees_scores(x, patient, cex = cex, file = f[2])
  revolver_plt_patient_trees(x, patient, max.phylogenies = max.phylogenies, 
                             cex = cex, file = f[3])

  ch = ceiling(sqrt(max.phylogenies))

  jamPDF(f[1:2], out.file = paste0(patient, '-data_scores.pdf'), layout = '2x1',crop.white = TRUE, page = 'a4')
  jamPDF(f[3], f[3], layout = paste0(ch,'x',ch), crop.white = TRUE,  page = 'a4')
  jamPDF(c(paste0(patient, '-data_scores.pdf'), f[3]), file, layout = '1x2', crop.white = TRUE, page = 'a4')

  # wat
  invisible(NULL)
}
