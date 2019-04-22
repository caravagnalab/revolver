#' @title Generate PDF report for a patient's fit
#'
#' @param x An object of class \code{"rev_cohort_fit"}
#' @param patient The patient to plot
#' @param cex Cex for graphics
#' @param file Output file, cannot be NA.
#'
#' @import crayon
#'
#' @examples
#' data(Breast.fit)
#' revolver_report_fit_patient(Breast.fit, 'PD9770')
#'
#' @export
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
