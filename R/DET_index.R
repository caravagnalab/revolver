#' Title
#'
#' @param x 
#' @param drivers 
#' @param min.occurrences 
#'
#' Extract the information transfer available for a patient.
#' 
#' @description 
#' 
#' This function computes the index of Divergent Evolutionary 
#' Trajectories for the fit of the input cohort.
#' 
#' @param x A \code{REVOLVER} cohort with fits.
#' @param drivers The list of drivers to compute the DET for.
#' @param min.occurrences The minimum number of occurrences for
#' a trajectory to be considered, zero by default. See also 
#' function \code{\link{plot_penalty}}.
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
#' DET_index(TRACERx_NEJM_2017_REVOLVER)
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