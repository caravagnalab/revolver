# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# These functions are a series of getters that allow to query a cohort
# object - i.e,, a REVOLVER S3 object of type "rev_cohort", or any
# other S3 object that inherits from that class. These functions has been
# introduced along with a new implementation of the internal REVOLVER
# structure which is based on tibbles and the tidy paradigm.
#
# G. Caravagna. April 2019
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

#' Extract CCF/binary data for a patient
#'
#' @param x a REVOLVER object of class "rev_cohort" or "rev_cohort_fit"
#' @param p a patient
#'
#' @return CCF/binary data for \code{p}
#' @export
#'
#' @examples
#' data(CRC.cohort)
#' CCF(CRC.cohort, 'adenoma_3')
Data = function(x, p)
{
  x$dataset[[p]]
  # if(p %in% names(x$data)) return(x$data[[p]])
  #
  # stop('Input does not have data.frame entries for', p)
}

#' Extract CCF/binary driver data for a patient
#'
#' @param x a REVOLVER object of class "rev_cohort" or "rev_cohort_fit"
#' @param p a patient
#'
#' @return CCF/binary drivers data for \code{p}
#' @export
#'
#' @examples
#' data(CRC.cohort)
#' Drivers(CRC.cohort, 'adenoma_3')
Drivers = function(x, p){
  Data(x, p) %>% filter(is.driver)
}
# Drivers = Vectorize(Drivers, vectorize.args = 'p', SIMPLIFY = FALSE)

#' Extract CCF/binary data for truncal mutations of a patient
#'
#' @param x a REVOLVER object of class "rev_cohort" or "rev_cohort_fit"
#' @param p a patient
#'
#' @return CCF/binary data for \code{p}'s truncal mutations
#' @export
#'
#' @examples
#' data(CRC.cohort)
#' Drivers(CRC.cohort, 'adenoma_3')
Truncal = function(x, p){
  Data(x, p) %>% filter(is.clonal)
}


#' Return the id of the clonal cluster for this patient
#'
#' @param x a REVOLVER object of class "rev_cohort"
#' @param p a patient
#'
#' @return The id of the clonal cluster for \code{p}
#'
#' @export
#'
#' @examples
#' data(CRC.cohort)
#' Clonal_cluster(CRC.cohort, 'adenoma_3')
Clonal_cluster = function(x, p)
{
  CCF_clusters(x, p) %>% filter(is.clonal) %>% pull(cluster)
}



#' Extract CCF/binary data for subclonal mutations of a patient
#'
#' @param x a REVOLVER object of class "rev_cohort" or "rev_cohort_fit"
#' @param p a patient
#'
#' @return CCF/binary data for \code{p}'s subclonal mutations
#' @export
#'
#' @examples
#' data(CRC.cohort)
#' Subclonal(CRC.cohort, 'adenoma_3')
Subclonal = function(x, p){
  Data(x, p) %>% filter(!is.clonal)
}

#' Extract CCF/binary data for a patient
#'
#' @param x a REVOLVER object of class "rev_cohort" or "rev_cohort_fit"
#' @param p a patient
#'
#' @return CCF/binary data for \code{p}
#' @export
#'
#' @examples
#' data(CRC.cohort)
#' CCF(CRC.cohort, 'adenoma_3')
CCF = function(x, p)
{
  nbps = Samples(x,p)

  Data(x,p) %>%
    select(
      id,
      variantID,
      is.driver,
      is.clonal,
      cluster,
      !!nbps
    )
}

#' Extract samples name for a patient
#'
#' @param x a REVOLVER object of class "rev_cohort" or "rev_cohort_fit"
#' @param p a patient
#'
#' @return CCF/binary data for \code{p}
#' @export
#'
#' @examples
#' data(CRC.cohort)
#' CCF(CRC.cohort, 'adenoma_3')
Samples = function(x, p)
{
  names(x$CCF.parser(Data(x,p) %>% filter(row_number() == 1) %>% pull(CCF)))
}


#' Extract cluster-level CCF/binary data for a patient
#'
#' @param x a REVOLVER object of class "rev_cohort" or "rev_cohort_fit"
#' @param p a patient
#'
#' @return CCF/binary data for \code{p}
#' @export
#'
#' @examples
#' data(CRC.cohort)
#' CCF(CRC.cohort, 'adenoma_3')
CCF_clusters = function(x, p)
{
  x$CCF[[p]]
}



#' Extract phylogenetic or mutation trees for a patient
#'
#' @param x a REVOLVER object of class "rev_cohort" or "rev_cohort_fit"
#' @param p a patient
#'
#' @return phylogenetic or mutation trees data for \code{p}
#' @export
#'
#' @examples
#' data(CRC.cohort)
#' Phylo(CRC.cohort, 'adenoma_3')
Phylo = function(x, p, rank = NULL)
{
  # if(is.null(x$phylogenies[[p]])) stop('There are no phylogenies for ', p)
  if(is.null(rank)) x$phylogenies[[p]]
  else x$phylogenies[[p]][[1]]
}

# Phylo = Vectorize(Phylo, vectorize.args = 'p', SIMPLIFY = FALSE)


#' Extract the information transfer for a patient's tree
#'
#' @param x A REVOLVER cohort
#' @param p A patient
#' @param rank The rank of the tree to extract
#' @param type Either `clones` or `drivers`.
#'
#' @return The information transfer for a patient's tree
#' @export
#'
#' @examples
#' data(CRC.cohort)
#' Phylo(CRC.cohort, 'adenoma_3')
ITransfer = function(x, p, rank = 1, type = 'drivers')
{
  tree = Phylo(x, p, rank)

  if(type == 'drivers') return(tree$transfer$drivers)
  if(type == 'clones') return(tree$transfer$clones)
}


#' Extract fitted model for a patient
#'
#' @param x a REVOLVER object of class "rev_cohort_fit"
#' @param p a patient
#'
#' @return fitted model for  \code{p}
#' @export
#'
#' @examples
#' data(CRC.cohort)
#'
#' Fit(revolver_fit(CRC.cohort), 'adenoma_3')
Fit = function(x, p)
{
  if(is.null(x$fit$phylogenies[[p]])) stop('There is no fit for ', p)

  x$fit$phylogenies[[p]]
}

#' Return summary stastics for the cohort's patients
#'
#' @description Returns the number of samples per patient, the number
#' of drivers, the number of clonal and subclonal mutations etc. The
#' function can be run on a subset of patients.
#'
#' @param x A REVOLVER cohort.
#' @param patients The patients for which the summaries are computed.
#'
#' @return A tibble with summary stastics.
#'
#' @export
#'
#' @examples
#' data(CRC.cohort)
#' Stats(CRC.cohort)
Stats = function(x, patients = x$patients)  {

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

#' Return summary stastics for the cohort's drivers
#'
#' @description Returns the number of clonal and subclonal occurrences
#' of the drivers in the cohort, and their percentage relative to the
#' cohort size. The function can be run on a subset of drivers.
#'
#' @param x a REVOLVER object of class "rev_cohort"
#' @param drivers The drivers for which the summaries are computed.
#'
#' @return A tibble with the driver stastics.
#'
#' @export
#'
#' @examples
#' data(CRC.cohort)
#' Stats_drivers(CRC.cohort)
Stats_drivers = function(x, drivers = x$variantIDs.driver) {

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

# ' Return summary stastics for the cohort's trees
#'
#' @description Returns the number of clonal and subclonal occurrences
#' of each driver in the cohort, and their percentage relative to the
#' cohort size. The function can be run on a subset of patients.
#'
#' @param x a REVOLVER object of class "rev_cohort"
#' @param patients The patients for which the summaries are computed.
#'
#' @return A tibble with the driver stastics.
#'
#' @export
#'
#' @examples
#' data(CRC.cohort)
#' Stats_drivers(CRC.cohort)
Stats_trees = function(x, patients = x$patients) {

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


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# These are non exported getters that help coding
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

has_patient_trees = function(x, p = NULL)
{
  if (is.null(p))
    return(!is.null(x$phylogenies))
  else
    return(!is.null(x$phylogenies[[p]]))
}



# obj_has_trees = function(x)
# {
#   if(is.null(x$phylogenies))
#     stop('Cannot proceed without having computed the trees -- are you calling the right function?')
# }
