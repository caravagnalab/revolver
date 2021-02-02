#' Compute CCF-based clone trees for a patient (phylogenies).
#'
#' @description Create clone tree phylogenies from the package
#' \code{ctree}, fitting the input data given to \code{REVOLVER}.
#' A set of patient ids can be given as input, by default all are
#' used. A parmeter controls if you want to overwrite the trees
#' of the patient, where they already computed. Please refer to the
#' parameters of the \code{ctrees} constructor rom package \code{ctree}
#' for the input parameters that you can use to tune the construction.
#'
#' @param x A \code{REVOLVER} cohort.
#' @param patients A set of patient ids in the cohort, for which the
#' phylogenies are created.
#'
#' @family Cohort creation
#'
#' @param x A \code{REVOLVER} cohort with now available
#' phylogeneis for the required \code{patients}.
#'
#' @export
#'
#' @import ctree
#'
#' @examples
#' # Data released in the 'evoverse.datasets'
#' data('TRACERx_NEJM_2017_REVOLVER', package = 'evoverse.datasets')
#'
#' print(TRACERx_NEJM_2017_REVOLVER)
#'
#' # We use the standard parameters with overwrite = TRUE
#' # otherwise the computation skips the patient
#' TRACERx_NEJM_2017_REVOLVER = compute_clone_trees(
#'    TRACERx_NEJM_2017_REVOLVER,
#'    patients = "CRUK0002",
#'    overwrite = TRUE)
#'
#' Phylo(TRACERx_NEJM_2017_REVOLVER, 'CRUK0002', rank = 1)
compute_clone_trees = function(
  x,
  patients = Stats_cohort(x) %>% pull(patientID),
  overwrite = FALSE,
  ...
)
{
  stop_not_revolver_object(x)
  lapply(patients, stop_invalid_patient, x = x)

  if(!has_patient_trees(x)) x$phylogenies = NULL

  cli::cli_h1("Constructing Clone Trees [ctree - https://caravagn.github.io/ctree/]")
  # pioTit("Constructing Clone Tree objects via ctree - https://caravagn.github.io/ctree/")
  
  # pioStr("Input patients.", suffix = '\n')
  # cat(paste0(patients, collapse = ', '), '\n')

  # Looping over all patients
  filos = easypar::run(
    FUN = function(patient)
    {
      cli::cli_h1("PatientID: {.field {patient}}")
      cat("\n")
      
      if(has_patient_trees(x, patient) & !overwrite)
      {
        message('Trees already available for ', patient, ', use overwrite = FALSE to force overwriting.')

        return(x$phylogenies[[patient]])
      }

      # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      # Run ctree constructor
      # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      obj = ctree::ctrees(
        CCF_clusters(x, patient),
        Drivers(x, patient),
        Samples(x, patient),
        patient = patient,
        ...
      )

      # Just show how many combinations we have
      comb = combination_of_information_transfer(x, patient)
      
      # pio::pioStr('\n Combinations of Information Transfer : ', comb, suffix = '\n')
      # cli::cl

      return(obj)
    },
    PARAMS = lapply(patients, list),
    parallel = FALSE,
    progress_bar = FALSE
  )
  
  # Check errors and notify if any of the cohort cannot be built
  if(length(filos) != length(patients)){
    cli::cli_alert_danger(
    "Errors computing trees. \\
    Check outputs and fix those, or remove the patients that raise it!")
    stop("Aborting.")
  }
  
  names(filos) = patients

  x$phylogenies = filos

  return(x)
}


