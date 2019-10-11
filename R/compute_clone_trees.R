#' Compute CCF-based phylogenies for a \code{REVOLVER} cohort.
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
#' @family Tree functions
#'
#' @return A modififed object of class \code{"rev_cohort"} with available
#' phylogeneis for \code{patients}.
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
#' # otherwise the tool returns retuns the trees
#' # already available
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
  
  pioTit("Constructing Clone Tree objects via ctree - https://caravagn.github.io/ctree/")
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
      obj = ctrees(
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
  names(filos) = patients
  
  x$phylogenies = filos

  return(x)
}


