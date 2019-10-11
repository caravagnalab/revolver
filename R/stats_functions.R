# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Functions for summary statistics
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

#' Summary statistics for the cohort's patients.
#'
#' @description Returns the number of samples per patient, the number
#' of drivers, the number of clonal and subclonal mutations etc. The
#' function can be run on a subset of patients.
#'
#' @param x A REVOLVER cohort.
#' @param patients The patients for which the summaries are computed.
#'
#' @family Summary statistics
#'
#' @return A tibble with summary stastics.
#'
#' @export
#'
#' @examples
#' # Data released in the 'evoverse.datasets'
#' data('TRACERx_NEJM_2017_REVOLVER', package = 'evoverse.datasets')
#'
#' # Get the stats for all patients
#' Stats(TRACERx_NEJM_2017_REVOLVER)
#'
#' # And subset the patients
#' Stats(TRACERx_NEJM_2017_REVOLVER, patients = c('CRUK0001', 'CRUK0002'))
Stats = function(x, patients = x$patients)
{
  stop_not_revolver_object(x)

  st = data.frame(
    patientID = patients,
    stringsAsFactors = FALSE
  )

  st$numBiopsies = sapply(st$patientID, function(w) length(Samples(x, w)))
  st$numMutations = sapply(st$patientID, function(w) nrow(Data(x, w)))

  st$numDriverMutations = sapply(st$patientID, function(w) nrow(Drivers(x, w)))
  st$numClonesWithDriver = sapply(st$patientID, function(w) length(unique(Drivers(x, w) %>% pull(cluster))))

  st$numTruncalMutations = sapply(st$patientID, function(w) nrow(Truncal(x, w)))
  st$numSubclonalMutations = sapply(st$patientID, function(w) nrow(Subclonal(x, w)))


  st %>% as_tibble()
}

#' Synonim for function \code{\link{Stats}}.
#'
#' @description Just a wrapper to the \code{\link{Stats}} function.
#'
#' @param ... Parameters forwarded to the \code{\link{Stats}}. function.
#'
#' @family Summary statistics
#'
#' @return See function \code{\link{Stats}}.
#'
#' @export
#'
#' @examples
#' # Data released in the 'evoverse.datasets'
#' data('TRACERx_NEJM_2017_REVOLVER', package = 'evoverse.datasets')
#'
#' # Get the stats for all patients
#' Stats_cohort(TRACERx_NEJM_2017_REVOLVER)
#'
#' # And subset the patients
#' Stats_cohort(TRACERx_NEJM_2017_REVOLVER, patients = c('CRUK0001', 'CRUK0002'))
Stats_cohort = function(...) { Stats(...) }

#' Summary statistics for the cohort's driver events.
#'
#' @description Returns the number of clonal and subclonal occurrences
#' of the drivers in the cohort, and their percentage relative to the
#' cohort size. The function can be run on a subset of drivers.
#'
#' @param x A REVOLVER cohort.
#' @param drivers The drivers for which the summaries are computed.
#'
#' @family Summary statistics
#'
#' @return A tibble with the driver stastics.
#'
#' @export
#'
#' @examples
#' # Data released in the 'evoverse.datasets'
#' data('TRACERx_NEJM_2017_REVOLVER', package = 'evoverse.datasets')
#'
#' # Get the stats for all patients
#' Stats_drivers(TRACERx_NEJM_2017_REVOLVER)
#'
#' # And subset the patients
#' Stats_drivers(TRACERx_NEJM_2017_REVOLVER, drivers = c('TP53', 'KRAS'))
Stats_drivers = function(x, drivers = x$variantIDs.driver) {

  stop_not_revolver_object(x)

  st = data.frame(
    variantID = drivers,
    row.names = drivers,
    stringsAsFactors = FALSE
  )

  clonal = sapply(
    x$patients,
    function(p) { CCF(x, p) %>% filter(is.driver, is.clonal, variantID %in% drivers) %>% pull(variantID) }
  )
  clonal = table(unlist(clonal))

  subclonal = sapply(
    x$patients,
    function(p) { CCF(x, p) %>% filter(is.driver, !is.clonal, variantID %in% drivers) %>% pull(variantID) }
  )
  subclonal = table(unlist(subclonal))

  st$numClonal = st$numSubclonal = 0

  st[names(clonal), 'numClonal'] = clonal
  st[names(subclonal), 'numSubclonal'] = subclonal

  st = st %>%
    as_tibble() %>%
    mutate(
      p_clonal = numClonal/x$n$patients,
      p_subclonal = numSubclonal/x$n$patients,
      N_tot = numClonal + numSubclonal,
      p_tot = N_tot / x$n$patient
    ) %>%
    select(
      variantID,
      numClonal, p_clonal,
      numSubclonal, p_subclonal,
      N_tot, p_tot
    ) %>%
    arrange (desc(numClonal), desc(numSubclonal))

  st %>% as_tibble()
}

#' Summary statistics for the cohort's trees.
#'
#' @description Returns the statistics about the trees that are available
#' in the cohort. The function can be run on a subset of patients.
#'
#' @param x A REVOLVER cohort.
#' @param patients The patients for which the summaries are computed.
#'
#' @family Summary statistics
#'
#' @return A tibble with the driver stastics.
#'
#' @export
#'
#' @examples
#' # Data released in the 'evoverse.datasets'
#' data('TRACERx_NEJM_2017_REVOLVER', package = 'evoverse.datasets')
#'
#' # Get the stats for all patients
#' Stats_trees(TRACERx_NEJM_2017_REVOLVER)
#'
#' # And subset the patients
#' Stats_trees(TRACERx_NEJM_2017_REVOLVER, patients = c('CRUK0001', 'CRUK0002'))
Stats_trees = function(x, patients = x$patients) {

  stop_not_revolver_object(x)

  if(
    !has_patient_trees(x) |
    !all(sapply(patients, has_patient_trees, x = x))
  )
    stop("There are no trees in this cohort object, or trees for some of the required patients are missing.
         Cannot compute this summary statistics, aborting.")

  st = data.frame(
    patientID = patients,
    row.names = patients,
    stringsAsFactors = FALSE
  )

  st$hasTrees = st$patientID %in% names(x$phylogenies)

  st$numTrees = sapply(st$patientID, function(p) length(x$phylogenies[[p]]))

  MSc = function(p) {
    if(p %in% names(x$phylogenies))
      return(max(sapply(seq_along(x$phylogenies[[p]]), function(t) x$phylogenies[[p]][[t]]$score)))
    NA
  }

  st$maxScore = sapply(st$patientID, MSc)

  mSc = function(p) {
    if(p %in% names(x$phylogenies))
      return(min(sapply(seq_along(x$phylogenies[[p]]), function(t) x$phylogenies[[p]][[t]]$score)))
    NA
  }

  st$minScore = sapply(st$patientID, mSc)

  st$combInfTransf = sapply(st$patientID, combination_of_information_transfer,  x = x)

  st %>% as_tibble()
}


#' Summary statistics for the cohort's fits.
#'
#' @description Returns a tibble that extends the result of
#' \code{\link{Stats_trees}} with information about the fit models.
#' Compared to summaries returns by other \code{Stats_*} functions,
#' the information from this one is precomputed.
#'
#' @param x A REVOLVER cohort where fits have been computed.
#' @param patients The patients for which the summaries are required.
#'
#' @family Summary statistics
#'
#' @return A tibble with the fits stastics.
#'
#' @export
#'
#' @examples
#' # Data released in the 'evoverse.datasets'
#' data('TRACERx_NEJM_2017_REVOLVER', package = 'evoverse.datasets')
#'
#' # Get the stats for all patients
#' Stats_fits(TRACERx_NEJM_2017_REVOLVER)
#'
#' # And subset the patients
#' Stats_fits(TRACERx_NEJM_2017_REVOLVER, patients = c('CRUK0001', 'CRUK0002'))
Stats_fits = function(x, patients = x$patients) {

  stop_not_revolver_object(x)

  if(
    !all(has_fits(x, patients))
  )
    stop("There are no fits in this cohort object, or fits for some of the required patients are missing.
         Cannot compute this summary statistics, aborting.")

  x$fit$fit_table %>%
    filter(patientID %in% patients)
}
