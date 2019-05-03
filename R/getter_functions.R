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
  stop_not_revolver_object(x)
  
  if(p %in% names(x$dataset)) return(x$dataset[[p]])
  stop(p, ' does does not have data in the REVOLVER cohort"', x$annotation, '"')
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
Drivers = function(x, p) 
{
  stop_not_revolver_object(x)
  
  Data(x, p) %>% filter(is.driver)
}
# Drivers = Vectorize(Drivers, vectorize.args = 'p', SIMPLIFY = FALSE)

#' Extract truncal mutations in a patient
#'
#' @param x A REVOLVER object of class "rev_cohort" or "rev_cohort_fit"
#' @param p A patient
#'
#' @return CCF/binary data for \code{p}'s truncal mutations
#' @export
#'
#' @examples
#' data(CRC.cohort)
#' Drivers(CRC.cohort, 'adenoma_3')
Truncal = function(x, p)
{
  stop_not_revolver_object(x)
  
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
#' data(TRACERx_cohort)
#' Clonal_cluster(TRACERx_cohort, 'CRUK0002')
Clonal_cluster = function(x, p)
{
  stop_not_revolver_object(x)
  
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
Subclonal = function(x, p)
{
  stop_not_revolver_object(x)
  
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
  stop_not_revolver_object(x)
  
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
#' @export Samples
#'
#' @examples
#' data(TRACERx_cohort)
#' Samples(TRACERx_cohort, 'CRUK0001')
Samples = function(x, p)
{
  stop_not_revolver_object(x)
  
  names(x$CCF_parser(Data(x,p) %>% filter(row_number() == 1) %>% pull(CCF)))
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
#' data(TRACERx_cohort)
#' CCF_clusters(TRACERx_cohort, 'CRUK0001')
CCF_clusters = function(x, p)
{
  stop_not_revolver_object(x)
  
  x$CCF[[p]]
}



#' Extract phylogenetic or mutation trees for a patient
#'
#' @param x a REVOLVER object of class "rev_cohort" or "rev_cohort_fit"
#' @param p a patient
#' @param data Either `trees` or `fit`.
#'
#' @return phylogenetic or mutation trees data for \code{p}
#' @export
#'
#' @examples
#' data(TRACERx_cohort)
#' 
#' # Get all the trees for a patient
#' Phylo(TRACERx_cohort, 'CRUK0002')
#' 
#' # Get a specific tree for a patient
#' Phylo(TRACERx_cohort, 'CRUK0002', rank = 2)
# 
#' # Get the fit tree for a patient
# Phylo(TRACERx_cohort, 'CRUK0002', data = 'fits')
Phylo = function(x, p, rank = NULL, data = 'trees')
{
  stop_not_revolver_object(x)
  
  if(!has_patient_trees(x, p, rank))
    stop('The requested tree does not exist, aborting.')

  if(data ==  'trees')
  {
    if(is.null(rank)) return(x$phylogenies[[p]])
    else {
      return(x$phylogenies[[p]][[rank]])
    }
  }
  
  if(data == 'fits')
    return(x$fit$phylogenies[[p]])
  
  stop("Unrecognized parameters, aborting")
}

# Phylo = Vectorize(Phylo, vectorize.args = 'p', SIMPLIFY = FALSE)


#' Extract the information transfer for a patient's tree
#'
#' @param x A REVOLVER cohort
#' @param p A patient
#' @param rank The rank of the tree to extract
#' @param type Either `clones` or `drivers`.
#' @param data Either `trees` or `fit`.
#'
#' @return The information transfer for a patient's tree
#' @export
#'
#' @examples
#' data(TRACERx_cohort)
#' 
#' # Get the transfer among drivers, top-ranking tree
#' Phylo(TRACERx_cohort, 'CRUK0002', rank = 1, type = 'drivers')
#' 
#' # Get the transfer among clones, top-ranking tree
#' Phylo(TRACERx_cohort, 'CRUK0002', rank = 1, type = 'clones')
#' 
#' # Get the transfer from the fit
#' Phylo(TRACERx_cohort, 'CRUK0002', rank = 1, type = 'clones', data = 'fits')
#' Phylo(TRACERx_cohort, 'CRUK0002', rank = 1, type = 'drivers', data = 'fits')
ITransfer = function(x, p, rank = 1, type = 'drivers', data = 'trees')
{
  stop_not_revolver_object(x)
  
  if(data ==  'trees')
  {
    tree = Phylo(x, p, rank)

    if(type == 'drivers') return(tree$transfer$drivers)
    if(type == 'clones') return(tree$transfer$clones)
  }
  
  if(data ==  'fits')
  {
    if(rank != 1) stop("If you ask for fit data rank must be 1, aborting.")
  
    tree = x$fit$phylogenies[[p]]
    
    if(type == 'drivers') return(tree$transfer$drivers)
    if(type == 'clones') return(tree$transfer$clones)
  }
  
  stop("Unrecognized parameters, aborting")
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
# Fit = function(x, p)
# {
#   if(is.null(x$fit$phylogenies[[p]])) stop('There is no fit for ', p)
# 
#   x$fit$phylogenies[[p]]
# }

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Clusters
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


#' Get clustering information for patients.
#' 
#' @description For a set of `patients` - all by default - 
#' extract the clustering information computed by REVOVLER.
#'
#' @param x A REVOLVEr cohort with fits and clusters.
#' @param patients Patients to use, all by default.
#'
#' @return A tibble with the required clustering assignments
#' 
#' @export
#'
#' @examples
#' TODO
Cluster = function(x, patients = x$patients)
{
  stop_not_revolver_object(x)
  
  if(!all(has_clusters(x, patients)))
    stop("There is no cluster information for what you're looking for,  aborting.")
  
  data.frame(
    patientID = patients,
    cluster = x$cluster$fits$labels[patients], 
    stringsAsFactors = FALSE) %>% 
    as_tibble() 
}


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Functions for summary statistics
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

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
#' data(TRACERx_cohort)
#' 
#' # Get the stats for all patients
#' Stats(TRACERx_cohort) 
#' 
#' # And subset the patients
#' Stats(TRACERx_cohort, patients = c('CRUK0001', 'CRUK0002')) 
Stats = function(x, patients = x$patients)  {

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

#' See \code{Stats} function.
#'
#' @description Just a wrapper to the \code{Stats} function.
#'
#' @param ... Parameters forwarded to the \code{Stats} function.
#'
#' @return A tibble with summary stastics.
#'
#' @export
#'
#' @examples
#' data(TRACERx_cohort)
#' 
#' # Get the stats for all patients
#' Stats_cohort(TRACERx_cohort) 
#' 
#' # And subset the patients
#' Stats_cohort(TRACERx_cohort, patients = c('CRUK0001', 'CRUK0002')) 
Stats_cohort = function(...) { Stats(...) }

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
#' data(TRACERx_cohort)
#' 
#' # Get the stats for all patients
#' Stats_drivers(TRACERx_cohort) 
#' 
#' # And subset the patients
#' Stats_drivers(TRACERx_cohort, patients = c('CRUK0001', 'CRUK0002')) 
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

#' Return summary stastics for the cohort's trees
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
#' data(TRACERx_cohort)
#' 
#' # Get the stats for all patients
#' Stats_trees(TRACERx_cohort) 
#' 
#' # And subset the patients
#' Stats_trees(TRACERx_cohort, patients = c('CRUK0001', 'CRUK0002')) 
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


#' Return summary stastics for the cohort's fits
#'
#' @description Returns a tibble that extends the result of 
#' \code{Stats_trees} with information about the fit models.
#' Compared to summaries returns by other \code{Stats_*} functions,
#' the information from this one is precomputed.
#'
#' @param x A REVOLVER cohort object with fits.
#' @param patients The patients for which the summaries are required.
#'
#' @return A tibble with the fits stastics.
#'
#' @export
#'
#' @examples
#' data(TRACERx_cohort)
#' 
#' # Get the stats for all patients
#' Stats_fits(TRACERx_cohort) 
#' 
#' # And subset the patients
#' Stats_fits(TRACERx_cohort, patients = c('CRUK0001', 'CRUK0002')) 
Stats_fits = function(x, patients = x$patients) {
  
  stop_not_revolver_object(x)
  
  if(
    !has_fits(x) |
    !all(sapply(patients, has_fits, x = x))
    ) 
    stop("There are no fits in this cohort object, or fits for some of the required patients are missing. 
         Cannot compute this summary statistics, aborting.")
  
  x$fit$fit_table %>%
    filter(patientID %in% patients)
  }

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# To test if the object has some consistency internally
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

stop_not_revolver_object = function(x)
{
  if(!inherits(x, 'rev_cohort'))
    stop("Input object is not a REVOLVER cohort, aborting.")
}

has_patient_trees = function(x, p = NULL, rank = NULL)
{
  if (is.null(p))
    return(!is.null(x$phylogenies))

  if (is.null(rank))
    return(p %in% names(x$phylogenies))
  else
    return(rank <= length(x$phylogenies[[p]]))
}

has_fits = function(x, p = NULL)
{
  if (is.null(p))
    return(!is.null(x$fit))
  
  return(!is.null(x$fit$phylogenies[[p]]))
}
has_fits = Vectorize(has_fits, vectorize.args = 'p', SIMPLIFY = TRUE)  

has_clusters= function(x, p = NULL)
{
  if (is.null(p))
    return(!is.null(x$cluster))
  
  return(p %in% names(x$cluster$fits$labels))
}
