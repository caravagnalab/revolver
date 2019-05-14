plot_drivers_map = function(x, cex = 1)
{
  # Get drivers, mirrored
  st = Stats_drivers(x) %>%
    arrange(desc(numClonal), desc(numSubclonal))

  order_drivers = st$variantID
  st$numSubclonal = - st$numSubclonal

  st = st %>%
    select(variantID, numClonal, numSubclonal) %>%
    rename(Clonal = numClonal, Subclonal = numSubclonal) %>%
    reshape2::melt(id = 'variantID')

  st$variantID = factor(st$variantID, levels = order_drivers)

  # Percentages
  N = length(x$patients)

  ggplot(st,
         aes(x = variantID,
             y = value,
             fill = variable)) +
    geom_bar(stat="identity", position="identity") +
    coord_flip(clip = 'off') +
    scale_fill_manual(values = c(`Clonal` = 'steelblue', `Subclonal` = 'darkorange')) +
    theme_minimal(base_size = 10 * cex) +
    theme(
      legend.position = 'bottom',
      legend.key.size = unit(3, 'mm')
    ) +
    labs(
      title = paste("Driver burden"),
      y = paste0('Occurrences (n = ', N, ' patients)'),
      x = 'Driver',
      subtitle = paste(x$annotation)
    ) +
    guides(
      fill = guide_legend("")
    )
}
