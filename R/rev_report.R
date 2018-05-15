
#' @title Generate PDF report for a patient's data and trees
#'
#' @param x An object of class \code{"rev_cohort"}
#' @param patient The patient to plot
#' @param cex Cex for graphics
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
                                   file = paste0('REVOLVER-report-patient-', patient, '.pdf'))
{
  obj_has_trees(x)
  pio::pioHdr(paste('REVOLVER Plot: report on data and trees for', patient), toPrint = NULL)

  if(is.na(file)) stop("Output file cannot be NA")

  max.phylogenies = 20

  revolver_plt_patient_data(x, patient, cex = cex, file = '1.pdf')
  revolver_plt_patient_trees_scores(x, patient, cex = cex, file = '2.pdf')
  revolver_plt_patient_trees(x, patient, max.phylogenies = max.phylogenies, cex = cex, file = '3.pdf')

  ch = round(sqrt(max.phylogenies))
  # xx =  jamPDF('3.pdf', out.file = '3.pdf', layout = paste0(ch,'x',ch), hide.output = TRUE, delete.original = TRUE)

  # xx =  jamPDF(c('1.pdf', '2.pdf'), out.file = '12.pdf', layout = '2x1', page = 'a4')
  # xx =  jamPDF(c('12.pdf', '3.pdf'), out.file = '12.pdf', layout = '2x1', page = 'a4')
  # xx =  jamPDF(c('1.pdf', '2.pdf', '3.pdf'), out.file = 'K.pdf', layout = '3x1', page = 'a4', hide.output = FALSE, delete.original = FALSE)
  xx =  jamPDF(c('1.pdf', '2.pdf', '3.pdf'), out.file = file, layout = '1x1', hide.output = TRUE, delete.original = T, crop.white = F)
  xx =  jamPDF(file, out.file = file, layout = '4x4', hide.output = TRUE, delete.original = T, crop.white = F)

  invisible(NULL)
}



#' @title Generate PDF report for a patient's fit
#'
#' @param x An object of class \code{"rev_cohort_fit"}
#' @param patient The patient to plot
#' @param cex Cex for graphics
#' @param file Output file, cannot be NA.
#'
#' @return nothing
#' @export
#' @import crayon
#'
#' @examples
#' data(Breast.fit)
#' revolver_report_fit_patient(Breast.fit, 'PD9770')
revolver_report_fit_patient =  function(x,
                                        patient,
                                        cex = 1,
                                        file = paste0('REVOLVER-report-fit-patient-', patient, '.pdf'))
{
  obj_has_trees(x)
  pio::pioHdr(paste('REVOLVER Plot: report on fit for', patient), toPrint = NULL)

  if(is.na(file)) stop("Output file cannot be NA")

  revolver_plt_fit_patient(x, patient, cex = cex, file = '1.pdf')
  revolver_plt_trajectories_patient(x, patient, cex = cex, file = '2.pdf')
  revolver_plt_itransfer_patient(x, patient, cex = cex, file = '3.pdf')

  xx =  jamPDF(c('1.pdf', '2.pdf', '3.pdf'), out.file = file, layout = '3x1')

  invisible(NULL)
}



# quartz()
# par(mfrow = c(1, 3))
# revolver_plt_fit_patient(Breast.fit, patient = 'PD9850')
# revolver_plt_trajectories_patient(Breast.fit, patient = 'PD9850')
# revolver_plt_itransfer_patient(Breast.fit, patient = 'PD9850')
