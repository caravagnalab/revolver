#' Check basic inconsistencies in a REVOLVER cohort.
#'
#' @description  Perform some basic diagnostic of a cohort object.
#' It will inform of patients without drivers and other
#' information that can be used to reshape the data before
#' fitting a model.
#'
#' @param x A \code{REVOLVER} cohort.
#' @param stopOnError Whether or not it should raise a stop on error.
#'
#' @return Nothing, all relevant information are print to screen
#' and the stop is raised only if \code{stopOnError = TRUE}.
#' @family Cohort creation
#' @export
#'
#' @examples
#' # Data released in the 'evoverse.datasets'
#' data('TRACERx_NEJM_2017_REVOLVER', package = 'evoverse.datasets')
#' 
#' revolver_check_cohort(TRACERx_NEJM_2017_REVOLVER)
#' 
#' print(TRACERx_NEJM_2017_REVOLVER) # calls this anyway
revolver_check_cohort = function(x, stopOnError = FALSE)
{
  # Duplicate IDs
  duplicate_varIDs = get_duplicates(x)
  
  
  if(nrow(duplicate_varIDs) > 0)
  {
    cli::boxx("ERROR - Duplicated driver variantIDs in the cohort should be removed!", 
              padding = 1, 
              background_col = "brown", 
              col = 'white', 
              float = 'center') %>% 
      cat()
    
    cat("\n")
    print(duplicate_varIDs)

    if(stopOnError) stop("Aborting.")
  }

  # Non-recurring drivers
  uncorr = Stats_drivers(x) %>%
    filter(N_tot <= 1)

  if(nrow(uncorr) > 0)
  {
    cli::boxx("WARNING - Driver variantIDs occuring only once could be removed.", 
              padding = 1, 
              background_col = "orange", 
              col = 'black', 
              float = 'center') %>% 
      cat()
    
    cat("\n")
    
    # cat(crayon::bgYellow('\n WARNING '),
    #     red("Some driver variantIDs occur only once and should therefore be removed. \n\n\tYou can use `revolver::Stats_drivers(x)` to retrieve them."), '\n')
    print(uncorr)

    if(stopOnError) stop("Aborting.")
  }

  # Patients with no drivers
  noDrv = Stats(x) %>%
    filter(numDriverMutations == 0)

  if(nrow(noDrv) > 0)
  {
    cli::boxx("WARNING - Patients without drivers could be removed.", 
              padding = 1, 
              background_col = "orange", 
              col = 'black', 
              float = 'center') %>% 
      cat()
    
    cat("\n")
    
    # cat(crayon::bgYellow('\n WARNING '),
    #     red("Some patients have no drivers and should therefore be removed. \n\n\tYou can use `revolver::Stats(x)` to retrieve them."), '\n')
    # print(noDrv)

    if(stopOnError) stop("Aborting.")
  }

  # Patients with 1 clone with drivers
  oneClDrv = Stats(x) %>%
    filter(numClonesWithDriver == 1)

  if(nrow(oneClDrv) > 0)
  {
    cli::boxx("WARNING - Some patients have only one clone with drivers; they will just be expanded.", 
              padding = 1, 
              background_col = "orange", 
              col = 'black', 
              float = 'center') %>% 
      cat()
    
    cat("\n")
    
    # cat(crayon::bgYellow('\n WARNING '),
    #     red("Some patients have only one clone with drivers, and therefore they will just be expanded."), '\n')
    print(oneClDrv)
  }

}

# Return duplicated IDs for drivers
get_duplicates = function(x)
{
  Reduce(
    bind_rows,
    lapply(
      x$patients,
      function(patient)
      {
        Data(x, patient) %>%
          filter(is.driver) %>%
          group_by(variantID) %>%
          summarise(occurrences = n()) %>%
          filter(occurrences >1) %>%
          arrange(desc(occurrences)) %>%
          ungroup() %>%
          mutate(patientID = patient)
      })
  )
}
