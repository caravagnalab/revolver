#' Compute binary mutation trees for a patient (phylogenies).
#'
#' @description Create mutation tree phylogenies from the package
#' \code{btree}, fitting the input data given to \code{REVOLVER}.
#' A set of patient ids can be given as input, by default all are
#' used. A parmeter controls if you want to overwrite the trees
#' of the patient, where they already computed. Please refer to the
#' parameters of the \code{btree} constructor from package \code{btree}
#' for the input parameters that you can use to tune the construction.
#' 
#' @note This function is not yet implemented
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
#' @import mtree
#'
#' @examples
#' print('TODO')
compute_mutation_trees = function(
  x,
  patients = Stats_cohort(x) %>% pull(patientID),
  overwrite = FALSE,
  ...
)
{
  stop_not_revolver_object(x)
  lapply(patients, stop_invalid_patient, x = x)

  if(!has_patient_trees(x)) x$phylogenies = NULL

  pioTit("Constructing Mutations Tree objects via mtree - https://caravagn.github.io/mtree/")
  pioStr("Input patients.", suffix = '\n')
  cat(paste0(patients, collapse = ', '), '\n')
  
  # Looping over all patients
  filos = easypar::run(
    FUN = function(patient)
    {
      if(has_patient_trees(x, patient) & !overwrite)
      {
        message('Trees already available for ', patient, ', use overwrite = FALSE to force overwriting.')

        return(x$phylogenies[[patient]])
      }

      # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      # Run ctree constructor
      # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      obj = mtree::mtrees(
        CCF_clusters(x, patient),
        Drivers(x, patient),
        Samples(x, patient),
        patient = patient,
        ...
      )

      # Just show how many combinations we have
      comb = combination_of_information_transfer(x, patient)
      pio::pioStr('\n Combinations of Information Transfer : ', comb, suffix = '\n')

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


