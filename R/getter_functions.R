### Getters

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
Phylo = function(x, p)
{
  if(is.null(x$phylogenies[[p]])) stop('There are no phylogenies for ', p)

  x$phylogenies[[p]]
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


# CCF = Vectorize(CCF, vectorize.args = 'p')
# Phylo = Vectorize(Phylo, vectorize.args = 'p')
# Fit = Vectorize(Fit, vectorize.args = 'p')
#
# TRACERx.cohort = fit
# TRACERx.cohort$fit = NULL
# class(TRACERx.cohort) = "rev_cohort"
# save(TRACERx.cohort, file = 'TRACERx.cohort.RData')
#
# TRACERx = TRACERx.cohort$dataset
# head(TRACERx)
#
# save(TRACERx, file = 'TRACERx.RData')
#

#
# CCF(x, p)
# Phylo(x, p)
# Fit(x, p)
# Fit(x, c(p,p))

# x= fit
# p = 'CRUK0001'
# load("/Users/gcaravagna/Documents/GitHub/test.revolver/TRACERx-release/Drivers_TabS23_in_>=2_patients/TRACERx.fit.RData")

Stats = function(x) {

  st = data.frame(
    patientID = x$patients,
    stringsAsFactors = FALSE
  )

  st$numBiopsies = sapply(st$patientID, function(w) length(Samples(x, w)))
  st$numMutations = sapply(st$patientID, function(w) nrow(Data(x, w)))
  st$numDriverMutations = sapply(st$patientID, function(w) nrow(Drivers(x, w)))
  st$numTruncalMutations = sapply(st$patientID, function(w) nrow(Truncal(x, w)))
  st$numSubclonalMutations = sapply(st$patientID, function(w) nrow(Subclonal(x, w)))


  st %>% as_tibble()
}

