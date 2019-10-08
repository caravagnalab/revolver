#' Plot the oncoprint for a patient.
#'
#' @param x A REVOLVER cohort object.
#' @param patient The id of a patient.
#' @param clusters_palette A palette function that should return the colour of
#' an arbitrary number of clusters.
#' @param ... Extra parameters, not used.
#'
#' @return A figure assembled with \code{ggpubr} which combines multiple
#' \code{ggplot} tile plots.
#' 
#' @export
#' 
#'
#' @examples
#' TODO
plot_patient_oncoprint = function(x,
                     patient,
                     clusters_palette = revolver:::distinct_palette_many,
                     ...)
{
  # Samples and CCF are all we use to make this plot
  samples = Samples(x, patient)
  values = Data(x, patient) %>%
    select(!!samples, id, cluster, is.clonal)
  
  values$sumCCF = apply(values[, samples, drop = FALSE],
                        1,
                        sum)
  
  # Cluster-specific ordering
  CCF_clusters = CCF_clusters(x, patient) %>%
    select(cluster,!!samples)
  
  CCF_clusters$sumCCF_cl = apply(CCF_clusters[, samples, drop = FALSE],
                                 1,
                                 function(w)
                                   2 ^ {
                                     sum(w > 0)
                                   } * sum(w))
  
  CCF_clusters = CCF_clusters %>% arrange(desc(sumCCF_cl))
  
  # Order according to a cobination of scores
  values = values %>%
    left_join(CCF_clusters %>% select(cluster, sumCCF_cl),
              by = 'cluster') %>%
    arrange(desc(is.clonal), desc(sumCCF_cl), desc(sumCCF))
  
  ordering_mutations = values$id
  
  # Melt everything and make the first plot
  CCF_values = values %>%
    select(!!samples, id) %>%
    reshape2::melt(id = 'id')
  CCF_values$id = factor(CCF_values$id, level = ordering_mutations)
  
  # Drivers to annotate
  driver_CCF = Drivers(x, patient) %>%
    select(variantID, !!samples[1], id) %>%
    reshape2::melt(id = c('id', 'variantID'))
  
  pl_CCF = ggplot(CCF_values,
                  aes(x = id, y = variable, fill = value)) +
    geom_tile(aes(height  = .9)) +
    scale_fill_distiller(palette = 'YlGnBu',
                         breaks = seq(0, 1, 0.2),
                         direction = 1) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(angle = 90, hjust = 0.5),
      legend.position = 'bottom',
      legend.key.size = unit(3, 'mm')
    ) +
    guides(fill = guide_colorbar("CCF",  barwidth = unit(5, 'cm'))) +
    labs(
      x = 'Mutation',
      y = 'Region') 
    # geom_text_repel(
    #   data = driver_CCF,
    #   aes(label = variantID),
    #   nudge_x = 0.5
    # )
  
  # Data for the second plot - clustering assignments
  Cluster_values = values %>%
    select(cluster, id) %>%
    reshape2::melt(id = 'id')
  Cluster_values$id = factor(Cluster_values$id, level = ordering_mutations)
  
  ncluster = length(unique(Cluster_values$value))
  
  # Second plot
  pl_Clusters = ggplot(Cluster_values,
                       aes(x = id, y = variable, fill = value)) +
    geom_tile(aes(height  = .9)) +
    scale_fill_manual(values = clusters_palette(ncluster)) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.text.y = element_text(angle = 90, hjust = 0.5),
      legend.position = 'top',
      legend.key.size = unit(3, 'mm')
    ) +
    guides(fill = guide_legend("Cluster", nrow = 1)) +
    labs(y = '')
  
  # Plot assembly
  require(ggpubr)
  
  figure = ggarrange(
    pl_Clusters,
    pl_CCF,
    heights = c(.3, 1),
    nrow = 2,
    ncol = 1
  )
  
  annotate_figure(
    figure,
    top = text_grob(
      label = paste0("Oncoprint for ", patient),
      hjust = 0)
  )
}
