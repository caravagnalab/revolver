
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

  revolver_plt_patient_data(x, patient, cex = cex, file = f[1])
  revolver_plt_patient_trees_scores(x, patient, cex = cex, file = f[2])
  revolver_plt_patient_trees(x, patient, max.phylogenies = max.phylogenies, cex = cex, file = f[3])

  ch = ceiling(sqrt(max.phylogenies))

  jamPDF(f[1:2], out.file = paste0(patient, '-data_scores.pdf'), layout = '2x1',crop.white = TRUE, page = 'a4')
  jamPDF(f[3], f[3], layout = paste0(ch,'x',ch), crop.white = TRUE,  page = 'a4')
  jamPDF(c(paste0(patient, '-data_scores.pdf'), f[3]), file, layout = '1x2', crop.white = TRUE, page = 'a4')

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

  f = paste(patient,
            c('fit.pdf', 'trajectories.pdf', 'itransfer.pdf'), sep = '-')

  revolver_plt_fit_patient(x, patient, cex = cex, file = f[1])
  revolver_plt_trajectories_patient(x, patient, cex = cex, file = f[2])
  revolver_plt_itransfer_patient(x, patient, cex = cex, file = f[3])

  jamPDF(f, out.file = file, layout = '3x1', page = 'a4', crop.white = TRUE)

  invisible(NULL)
}



# quartz()
# par(mfrow = c(1, 3))
# revolver_plt_fit_patient(Breast.fit, patient = 'PD9850')
# revolver_plt_trajectories_patient(Breast.fit, patient = 'PD9850')
# revolver_plt_itransfer_patient(Breast.fit, patient = 'PD9850')
