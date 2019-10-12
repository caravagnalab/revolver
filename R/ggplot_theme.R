my_ggplot_theme = function(cex = 1)
{
  theme_light(base_size = 10 * cex) +
    theme(
      legend.position = "bottom",
      legend.key.size = unit(.3 * cex, "cm"),
      panel.background = element_rect(fill = 'white')
    )
}
