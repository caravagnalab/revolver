
#' Plot the clusters data for a patient.
#' 
#' @description This function creates a \code{ggplot}-style
#' heatmap of the mean CCF of each clone in a patient's data.
#' The heatmap is annotated for the drivers status of each
#' clone (with/ without driver). The CCF values are used to
#' colour the heatmap (`NA` values are in light gray).
#'
#' @param x A REVOLVER cohort object.
#' @param patient A patient id.
#' @param cex Cex of the plot.
#' @param ... Extra parameters, not used.
#'
#' @return A \code{ggplot} plot.
#' @export
#'
#' @examples
#' TODO
plot_data_clusters = function(x, patient, cex = 1, ...)
{
  p = patient
  
  sm = Samples(x, p)
  cl = CCF_clusters(x, p)

  # Values CCF
  cl_tab = cl %>%
    select(!!sm, cluster) %>%
    reshape2::melt(id = 'cluster') %>%
    rename(
      region = variable,
      CCF = value) %>%
    as_tibble() %>%
    mutate(CCF = ifelse(CCF == 0, NA, CCF))
  
  # Annotations
  cl_tab_anno = cl %>%
     select(cluster, nMuts, is.driver, is.clonal)
  
  # Cluster ordering by sum of CCF
  cluster_ordering = cl_tab %>%
    group_by(cluster) %>%
    summarise(tot = sum(CCF, na.rm = TRUE)) %>%
    arrange(tot) %>%
    pull(cluster)
  
  # Combined
  cl_tab = cl_tab %>%
    left_join(cl_tab_anno %>% select(cluster, is.driver), by = 'cluster') %>%
    mutate(is.driver = ifelse(is.na(CCF), FALSE, is.driver))
  
  cl_CCF = CCF(x, p) %>%
    group_by(cluster) %>%
    filter(is.driver) %>%
    summarise(
      nDrivers = n(),
      label = paste(variantID, collapse = ', ')
    )
  
  cl_tab_anno = cl_tab_anno %>%
    left_join(cl_CCF, by = 'cluster') %>%
    mutate(
      nDrivers = ifelse(is.na(nDrivers), '0', nDrivers)
    )
  
  # Factor to sort
  cl_tab$cluster = factor(cl_tab$cluster, levels = cluster_ordering)
  # cl_CCF$cluster = factor(cl_CCF$cluster, levels = cluster_ordering)
  
  # Drivers number
  # mD = max(cl_CCF$nDrivers)
  # cl_tab_anno$nDrivers = factor(cl_tab_anno$nDrivers, levels = paste(0:mD))
  
  ggplot(
    cl_tab,
    aes(
      x = region,
      y = cluster,
      z = CCF,
      fill = CCF,
      color = is.driver)
  ) +
    geom_tile(aes(width = .8, height = .8), size = 1 * cex) +
    geom_text(aes(label = CCF), color = 'orange') +
    scale_fill_distiller(palette = 'Blues', na.value = 'gainsboro', direction = 1) +
    scale_color_manual(
      values = c(`TRUE` = 'red', `FALSE` = NA)
    ) +
    theme_minimal(base_size = 10 * cex) +
    theme(
      legend.position = 'bottom',
      legend.key.size = unit(3, 'mm')
    ) +
    labs(
      title = paste("Clones for", p),
      y = 'Cluster',
      x = 'Region',
      subtitle = paste("Clonal cluster = ", cl$cluster[cl$is.clonal])
    ) +
    guides(
      fill = guide_colourbar(barwidth = 6),
      color = guide_legend("With driver")
      )
}

