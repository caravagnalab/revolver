#' @title Check basic inconsistencies in a REVOLVER cohort.
#'
#' @details Perform some basic diagnostic of a cohort object.
#' It will inform of patients without drivers and other
#' information that can be used to reshape the data before
#' fitting a model.
#'
#' @param x An object of class \code{"rev_cohort"}
#' @param stopOnError Whether or not it should raise a stop on error.
#'
#' @return Nothing, all relevant information are print to screen
#' and the stop is raised only if \code{stopOnError = TRUE}.
#'
#' @import crayon
#'
#' @export
#'
#' @examples
#' data(CRC.cohort)
#' revolver_check_cohort(CRC.cohort)
#' print(CRC.cohort) # calls this anyway
revolver_check_cohort = function(x, stopOnError = FALSE)
{
  # Duplicate IDs
  duplicate_varIDs = get_duplicates(x)

  if(nrow(duplicate_varIDs) > 0)
  {
    cat(bgRed('\n ERROR '),
        red("Duplicated driver variantIDs in the cohort that should be removed. \n\n\tYou can use `revolver:::get_duplicates(x)` to retrieve them."), '\n')
    print(duplicate_varIDs)

    if(stopOnError) stop("Aborting.")
  }

  # Non-recurring drivers
  uncorr = Stats_drivers(x) %>%
    filter(N_tot <= 1)

  if(nrow(uncorr) > 0)
  {
    cat(bgRed('\n ERROR '),
        red("Some driver variantIDs occur only once and should therefore be removed. \n\n\tYou can use `revolver::Stats_drivers(x)` to retrieve them."), '\n')
    print(uncorr)

    if(stopOnError) stop("Aborting.")
  }

  # Patients with no drivers
  noDrv = Stats(x) %>%
    filter(numDriverMutations == 0)

  if(nrow(noDrv) > 0)
  {
    cat(bgRed('\n ERROR '),
        red("Some patients have no drivers and should therefore be removed. \n\n\tYou can use `revolver::Stats(x)` to retrieve them."), '\n')
    print(noDrv)

    if(stopOnError) stop("Aborting.")
  }

  # Patients with 1 clone with drivers
  oneClDrv = Stats(x) %>%
    filter(numClonesWithDriver == 1)

  if(nrow(oneClDrv) > 0)
  {
    cat(bgRed('\n WARNING '),
        red("Some patients have only one clone with drivers, and therefore they will just be expanded."), '\n')
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
