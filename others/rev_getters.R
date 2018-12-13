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
  x$dataset %>% filter(patientID %in% p)
  # if(p %in% names(x$data)) return(x$data[[p]])
  #
  # stop('Input does not have data.frame entries for', p)
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
  D = Data(x, p)

  foo = function(w)
  {
    v = revolver:::CCF.parser(w)
    tv = as_tibble(v)
    tv$sample = names(v)

    tv %>% spread(sample, value)
  }

  # CCF per row
  CCF = D %>%
    group_by(id) %>%
    do(foo(.$CCF)) %>%
    ungroup() %>%
    select(-id)

  CCF %>% mutate(cluster = D$cluster, id = D$id)
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
  colnames(CCF(x,p) %>% select(-cluster))
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
  pCCF = CCF(x,p)
  samples = colnames(pCCF)

  D = Data(x, p)

  nMuts = D %>% group_by(cluster) %>% summarise(nMuts = n())

  values = D %>% reshape2::melt(id = 'cluster') %>%
    as_tibble() %>%
    filter(variable %in% !!samples) %>%
    group_by(cluster) %>%
    summarise(median = median(as.numeric(value)))

  full_join(pCCF, nMuts, values, by = 'cluster') %>%
    select(cluster, nMuts, !!samples)
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
  if(is.null(x$phylogenies[[p]])) stop('Input does not have phylogenies for', p)

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
  if(is.null(x$fit$phylogenies[[p]])) stop('Input does not have a fit for', p)

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
