plot_cohort_drivers_patients = function(x)
{
  occurrences = lapply(
    x$patients, 
    Drivers, x = x)
  
  occurrences = Reduce(bind_rows, occurrences)
  
  occurrences = occurrences %>%
    select(patientID, variantID) 
  
  ggplot(
    occurrences,
    aes(
      y = variantID,
      x = patientID
    )
  ) +
    geom_tile() +
    theme(
      axis.text.x = element_text(angle = 90)
    )
  
}
