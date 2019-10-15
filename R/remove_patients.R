#' Remove patients from the cohort.
#'
#' @description Each patient is identied through its id (\code{patientID});
#' with this function, you can remove patients, which consists in removing  data
#' and trees. If you have fit the models or clustered the cohort,
#' you must re-run the analyses after this modification; for this reason,
#' any previous result from those analyses is cancelled from the returned object.
#'
#' Notice also that some drivers might be removed by this function.
#'
#' @param x A REVOLVER cohort.
#' @param patientID Id of the patient to remove. It can be a vector.
#'
#' @return A modified cohort without the required patients.
#' @export
#'
#' @examples
#' # Data released in the 'evoverse.datasets'
#' data('TRACERx_NEJM_2017_REVOLVER', package = 'evoverse.datasets')
#'
#' print(TRACERx_NEJM_2017_REVOLVER)
#'
#' new_cohort = remove_patient(TRACERx_NEJM_2017_REVOLVER, "CRUK0001")
#' print(new_cohort)
#' 
#' new_cohort = remove_patient(TRACERx_NEJM_2017_REVOLVER, c("CRUK0002", "CRUK00024"))
#' print(new_cohort)
remove_patients = function(x, 
                           patientID,
                           check = TRUE)
{
  new_patients = setdiff(x$patients, patientID)

  # Update counts and remove non-used patients
  x$patients = new_patients

  x$n$patients = length(x$patients)
  x$dataset = x$dataset[new_patients]
  x$CCF = x$CCF[new_patients]
  x$phylogenies = x$phylogenies[new_patients]

  # Because we removed drivers cancelling out this patient, we want to
  # check that there is no driver now occurring in 1 patient. If we
  # find it, we cancel if
  drv_to_cancel = Stats_drivers(x) %>% filter(N_tot <= 1) %>% pull(variantID)

  if(length(drv_to_cancel) > 0) {
    message(paste(drv_to_cancel, collapse = ', '), ": driver events that now are found in only one patient, will be now removed ...")
  }

   x = remove_driver(x, drv_to_cancel, check = FALSE)

  # if has clusters, force to recompute
  if(has_clusters(x)) {
    message("- The cohort has clusters which will be cancelled, please re-compute the clusters ...")
    x$cluster = NULL
  }

  # Check the cohort
  revolver_check_cohort(x)

  return(x)
}
