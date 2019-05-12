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