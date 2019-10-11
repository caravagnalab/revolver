#' Compute the index of Divergent Evolutionary Trajectories.
#' 
#' @description
#' 
#' The index of Divergent Evolutionary Trajectories is a 
#' measure derived from Shannon's entropy to determine,
#' for any driver event \code{X}, how heterogeneous are
#' the trajectories that lead to \code{X}.
#' 
#' This is simply based on counting the number of edges
#' \code{Y -> X}, for every \code{Y} in the data, and
#' using the distribution of the observed frequencies to
#' determine an entropy-derived measure (analog to 
#' species heterogeneity). This means that for values 
#' larger than 0, the model observers heterogeneous
#' trajectories in the data.
#' 
#' To compute the DET the cohort must have fits available.
#'
#' @param x A \code{REVOLVER} cohort object with fits available.
#' @param drivers The list of drivers to compute the DET for.
#' @param min.occurrences The minimum number of occurrences for
#' a trajectory to be considered, zero by default. 
#' 
#' @seealso Function \code{\link{plot_penalty}} plots another
#' measure derived from the same information used to compute the DET.
#' 
#' @family Summary statistics 
#' 
#' @return The DET index for the input cohort.
#' @export
#'
#' @examples
#' # Data released in the 'evoverse.datasets'
#' data('TRACERx_NEJM_2017_REVOLVER', package = 'evoverse.datasets')
#'  
#' # Get the DET with all cohort
#' DET_index(TRACERx_NEJM_2017_REVOLVER)
#' 
#' # Look specifically for TP53 - the DET suggests
#' # heterogeneous trajectories.
#' DET_index(TRACERx_NEJM_2017_REVOLVER, drivers = 'TP53')
#' 
#' # Look specifically for drivers in at least 5 patients
#' DET_index(TRACERx_NEJM_2017_REVOLVER, min.occurrences = 5)
DET_index = function(x, 
                     drivers = x$variantIDs.driver,
                     min.occurrences = 0
)
{
  # Subset E to make computations
  E = x$fit$penalty %>%
    filter(to %in% drivers, count >= min.occurrences)
  
  # Species counts 
  Species = E %>% 
    group_by(to) %>%
    summarise(N = n())
  
  Diversity = E %>% 
    group_by(to) %>%
    summarise(diversity = vegan::diversity(count))
  
  Diversity %>% 
    left_join(Species, by = 'to') %>%
    mutate(
      DET_index = diversity/log(N),
      DET_index = ifelse(is.na(DET_index), 0, DET_index)
    ) %>%
    rename(driver = to) %>%
    arrange(DET_index)
}