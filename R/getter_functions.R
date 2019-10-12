# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# These functions are a series of getters that allow to query a cohort
# object - i.e,, a REVOLVER S3 object of type "rev_cohort", or any
# other S3 object that inherits from that class. These functions has been
# introduced along with a new implementation of the internal REVOLVER
# structure which is based on tibbles and the tidy paradigm.
#
# G. Caravagna. October 2019
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


#' Extract all data for a patient.
#'
#' @description
#'
#' From a cohort object, this function extracts the tibble of
#' all data available for a patient (either CCF or binary).
#'
#' @param x A \code{REVOLVER} cohort.
#' @param p The id of a patient in the cohort.
#'
#' @family Getters
#'
#' @return Drivers data for a custom patient.
#' @export
#'
#' @examples
#' # Data released in the 'evoverse.datasets'
#' data('TRACERx_NEJM_2017_REVOLVER', package = 'evoverse.datasets')
#'
#' Data(TRACERx_NEJM_2017_REVOLVER, 'CRUK0002')
#'
#' Data(TRACERx_NEJM_2017_REVOLVER, 'CRUK0008')
Data = function(x, p)
{
  stop_not_revolver_object(x)

  if(p %in% names(x$dataset)) return(x$dataset[[p]])
  stop(p, ' does not have data in the REVOLVER cohort "', x$annotation, '"')
}

#' Extract driver data for a patient.
#'
#' @description
#'
#' From a cohort object, this function extracts the tibble of
#' annotated driver events for a patient (either CCF or binary).
#'
#' @param x A \code{REVOLVER} cohort.
#' @param p The id of a patient in the cohort.
#'
#' @family Getters
#'
#' @return Drivers data for a custom patient.
#' @export
#'
#' @examples
#' # Data released in the 'evoverse.datasets'
#' data('TRACERx_NEJM_2017_REVOLVER', package = 'evoverse.datasets')
#'
#' Drivers(TRACERx_NEJM_2017_REVOLVER, 'CRUK0002')
#'
#' Drivers(TRACERx_NEJM_2017_REVOLVER, 'CRUK0008')
Drivers = function(x, p)
{
  stop_not_revolver_object(x)

  Data(x, p) %>% filter(is.driver)
}
# Drivers = Vectorize(Drivers, vectorize.args = 'p', SIMPLIFY = FALSE)

#' Extract clonal (i.e., truncal) data for a patient.
#'
#' @description
#'
#' From a cohort object, this function extracts the tibble of
#' truncal events for a patient (either CCF or binary).
#'
#' @param x A \code{REVOLVER} cohort.
#' @param p The id of a patient in the cohort.
#'
#' @family Getters
#'
#' @return Truncal data for a custom patient.
#' @export
#'
#' @examples
#' # Data released in the 'evoverse.datasets'
#' data('TRACERx_NEJM_2017_REVOLVER', package = 'evoverse.datasets')
#'
#' Truncal(TRACERx_NEJM_2017_REVOLVER, 'CRUK0002')
#'
#' Truncal(TRACERx_NEJM_2017_REVOLVER, 'CRUK0008')
Truncal = function(x, p)
{
  stop_not_revolver_object(x)

  Data(x, p) %>% filter(is.clonal)
}

#' Extract the clonal cluster id for a patient.
#'
#' @description
#'
#' From a cohort object, this function return the
#' id of the clonal cluster for a patient.
#'
#' @param x A \code{REVOLVER} cohort.
#' @param p The id of a patient in the cohort.
#'
#' @family Getters
#'
#' @return The id of the clonal cluster for the patient.
#' @export
#'
#' @examples
#' # Data released in the 'evoverse.datasets'
#' data('TRACERx_NEJM_2017_REVOLVER', package = 'evoverse.datasets')
#'
#' Clonal_cluster(TRACERx_NEJM_2017_REVOLVER, 'CRUK0002')
#'
#' Clonal_cluster(TRACERx_NEJM_2017_REVOLVER, 'CRUK0008')
Clonal_cluster = function(x, p)
{
  stop_not_revolver_object(x)

  CCF_clusters(x, p) %>% filter(is.clonal) %>% pull(cluster)
}

#' Extract subclonal data for a patient.
#'
#' @description
#'
#' From a cohort object, this function extracts the tibble of
#' subclonal events for a patient (either CCF or binary).
#'
#' @param x A \code{REVOLVER} cohort.
#' @param p The id of a patient in the cohort.
#'
#' @family Getters
#'
#' @return Subclonal data for a custom patient.
#' @export
#'
#' @examples
#' # Data released in the 'evoverse.datasets'
#' data('TRACERx_NEJM_2017_REVOLVER', package = 'evoverse.datasets')
#'
#' Subclonal(TRACERx_NEJM_2017_REVOLVER, 'CRUK0002')
#'
#' Subclonal(TRACERx_NEJM_2017_REVOLVER, 'CRUK0008')
Subclonal = function(x, p)
{
  stop_not_revolver_object(x)

  Data(x, p) %>% filter(!is.clonal)
}

#' Extract all CCF or binary data for a patient.
#'
#' @description
#'
#' From a cohort object, this function extracts the tibble of
#' all CCF or binary data for a patient. This is just a faster
#' way to subset a generic call to \code{\link{Data}}.
#'
#' @param x A \code{REVOLVER} cohort.
#' @param p The id of a patient in the cohort.
#'
#' @family Getters
#'
#' @return Data for a custom patient.
#' @export
#'
#' @examples
#' # Data released in the 'evoverse.datasets'
#' data('TRACERx_NEJM_2017_REVOLVER', package = 'evoverse.datasets')
#'
#' CCF(TRACERx_NEJM_2017_REVOLVER, 'CRUK0002')
#'
#' CCF(TRACERx_NEJM_2017_REVOLVER, 'CRUK0008')
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

#' Extract all samples ids for a patient.
#'
#' @description
#'
#' From a cohort object, this function extracts a vector of
#' the sample ids available for a patient.
#'
#' @param x A \code{REVOLVER} cohort.
#' @param p The id of a patient in the cohort.
#'
#' @family Getters
#'
#' @return Sample names for a custom patient.
#' @export
#'
#' @examples
#' # Data released in the 'evoverse.datasets'
#' data('TRACERx_NEJM_2017_REVOLVER', package = 'evoverse.datasets')
#'
#' Samples(TRACERx_NEJM_2017_REVOLVER, 'CRUK0002')
#'
#' Samples(TRACERx_NEJM_2017_REVOLVER, 'CRUK0008')
Samples = function(x, p)
{
  stop_not_revolver_object(x)

  names(x$CCF_parser(Data(x,p) %>% filter(row_number() == 1) %>% pull(CCF)))
}

#' Extract all clusters in the data of a patient.
#'
#' @description
#'
#' From a cohort object, this function extracts a tibble of
#' all clusters in the data of a patient.
#'
#' @param x A \code{REVOLVER} cohort.
#' @param p The id of a patient in the cohort.
#'
#' @family Getters
#'
#' @return CCF or binary clusters data for a custom patient.
#' @export
#'
#' @examples
#' # Data released in the 'evoverse.datasets'
#' data('TRACERx_NEJM_2017_REVOLVER', package = 'evoverse.datasets')
#'
#' CCF_clusters(TRACERx_NEJM_2017_REVOLVER, 'CRUK0002')
#'
#' CCF_clusters(TRACERx_NEJM_2017_REVOLVER, 'CRUK0008')
CCF_clusters = function(x, p)
{
  stop_not_revolver_object(x)

  x$CCF[[p]]
}

#' Extract the trees available for a patient.
#'
#' @description
#'
#' From a cohort object, this function extracts the tree
#' objects that are available for a patient. Parameters
#' can be set to retrieve a particular tree, or the fit
#' tree which is however availble only after fitting the
#' cohort data. This function can either return a
#' phylogenetic clone tree (R object \code{ctree}), or mutation
#' trees  (R object \code{btree}).
#'
#' @param x A \code{REVOLVER} cohort.
#' @param p The id of a patient in the cohort.
#' @param rank The rank of the tree to extract.
#' @param data Either `trees` or `fit`, which requires to have
#' already computed the fit of the input cohort.
#'
#' @family Getters
#'
#' @return Phylogenetic or mutation trees data for the patient.
#'
#' @export
#'
#' @examples
#' # Data released in the 'evoverse.datasets'
#' data('TRACERx_NEJM_2017_REVOLVER', package = 'evoverse.datasets')
#'
#' # Get all the trees for a patient
#' Phylo(TRACERx_NEJM_2017_REVOLVER, 'CRUK0002')
#'
#' # Get a specific tree for a patient
#' Phylo(TRACERx_NEJM_2017_REVOLVER, 'CRUK0002', rank = 2)
#
#' # Get the fit tree for a patient
# Phylo(TRACERx_NEJM_2017_REVOLVER, 'CRUK0002', data = 'fits')
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


#' Extract the information transfer available for a patient.
#'
#' @description
#'
#' From a cohort object, this function extracts the tree
#' objects that are available for a patient. Parameters
#' can be set to retrieve a particular tree, or the fit
#' tree which is however availble only after fitting the
#' cohort data. This function can either return a
#' phylogenetic clone tree (R object \code{ctree}), or mutation
#' trees  (R object \code{btree}). Besides, it can return
#' the information transfer as the ordering of the annotated
#' drivers, or in terms o the ordering of the clones these
#' map to.
#'
#' @param x A \code{REVOLVER} cohort.
#' @param p The id of a patient in the cohort.
#' @param rank The rank of the tree to extract.
#' @param data Either `trees` or `fit`, which requires to have
#' already computed the fit of the input cohort.
#' @param type Either `drivers` or `clones`.
#'
#' @family Getters
#'
#' @return The information transfer for the tree of a patient
#' @export
#'
#' @examples
#' # Data released in the 'evoverse.datasets'
#' data('TRACERx_NEJM_2017_REVOLVER', package = 'evoverse.datasets')
#'
#' # Get the transfer among drivers, top-ranking tree
#' ITransfer(TRACERx_NEJM_2017_REVOLVER, 'CRUK0002', rank = 1, type = 'drivers')
#'
#' # Get the transfer among clones, top-ranking tree
#' ITransfer(TRACERx_NEJM_2017_REVOLVER, 'CRUK0002', rank = 1, type = 'clones')
#'
#' # Get the transfer from the fit
#' ITransfer(TRACERx_NEJM_2017_REVOLVER, 'CRUK0002', rank = 1, type = 'clones', data = 'fits')
#' ITransfer(TRACERx_NEJM_2017_REVOLVER, 'CRUK0002', rank = 1, type = 'drivers', data = 'fits')
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


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Clusters
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


#' Get REVOLVER clusters information for this cohort
#'
#' @description
#'
#' For a set of `patients` - all by default - extract the clustering information
#' computed by REVOVLER. The result is just a tibble with a column representing
#' the patient id, and a column with the cluster label.
#'
#' @param x A REVOLVEr cohort with fits and clusters.
#' @param patients Patients to use, all by default.
#'
#' @return A tibble with the required clustering assignments
#'
#' @export
#'
#' @examples
#' # Data released in the 'evoverse.datasets'
#' data('TRACERx_NEJM_2017_REVOLVER', package = 'evoverse.datasets')
#'
#' Cluster(TRACERx_NEJM_2017_REVOLVER)
#'
#' #' Cluster(TRACERx_NEJM_2017_REVOLVER, patients = c("CRUK0001", "CRUK0002"))
Cluster = function(x, patients = x$patients)
{
  stop_not_revolver_object(x)
  obj_has_fit(x)
  obj_has_clusters(x)

  data.frame(
    patientID = patients,
    cluster = x$cluster$fits$labels[patients],
    stringsAsFactors = FALSE) %>%
    as_tibble()
}


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# To test if the object has some consistency internally
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

stop_not_revolver_object = function(x)
{
  pass = any(sapply(c('rev_cohort', 'rev_cohort_fit'), inherits, x = x))

  if(!pass)
    stop("Input object is not a REVOLVER cohort, aborting.")
}

stop_invalid_patient = function(x, p)
{
  pass = p %in% x$patients

  if(!pass)
    stop("Patient", p, "does not exist in this REVOLVER cohort, aborting.")
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
