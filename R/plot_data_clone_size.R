#' Plot the clone size data for a patient.
#' 
#' @description This function creates a \code{ggplot}
#' barplot of the clone size values foe each clone in a patient's data.
#' The size of a clone is defined as the number of mutations assigned 
#' to it. The barplot is annotated to report wether a subclone with a driver
#' is significantly larger than the expected size for a subclone without
#' driver. To carry out this test subclones without drivers are used to
#' estimate the parameters of a univariate Gaussian distribution (mean
#' and standard deviation), the p-value is then computed from the fit
#' distribution through the `pnorm` function. The confidence level for
#' the test can be passed as parameter.
#'
#' @param x A REVOLVER cohort object.
#' @param patient A patient id.
#' @param alpha_level Alpha level for the test, default is 0.05.
#' @param cex Cex of the plot.
#' @param ... Extra parameters, not used.
#'
#' @return A \code{ggplot} plot.
#' 
#' @export
#'
#' @examples
#' TODO
plot_data_clone_size = function(x, patient, alpha_level = 0.05, cex = 1, ...)
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
  cl_tab_anno$cluster = factor(cl_tab_anno$cluster, levels = cluster_ordering)
  cl_CCF$cluster = factor(cl_CCF$cluster, levels = cluster_ordering)
  
  # Drivers number
  mD = max(cl_CCF$nDrivers)
  cl_tab_anno$nDrivers = factor(cl_tab_anno$nDrivers, levels = paste(0:mD))
  
  # Carry out the subclone size test
  pvals = pval_subcsz(x, p)
  
  cl_tab_anno = cl_tab_anno %>%
    mutate(cluster = paste(cluster)) %>%
    left_join(pvals, by = 'cluster') %>%
    mutate(
      significant = ifelse(pvalue < !!alpha_level, "YES", "NO"),
    )
  cl_tab_anno$cluster = factor(cl_tab_anno$cluster, levels = cluster_ordering)
  
  ggplot(
    cl_tab_anno,
    aes(
      x = cluster,
      y = nMuts,
      fill = nDrivers,
      color = significant
    )
  ) +
    geom_bar(stat = 'identity') +
    coord_flip(clip = 'off') +
    scale_fill_brewer(palette = 'Purples') +
    scale_color_manual(values = c(`YES` = 'forestgreen', `NO` = 'red')) +
    geom_text(aes(label = label), color = 'black', size = 2.5  * cex, hjust = 0) +
    theme_minimal(base_size = 10  * cex) +
    theme(
      legend.position = 'bottom',
      legend.key.size = unit(3, 'mm')
    ) +
    labs(
      title = paste0("Clones size for ", patient),
      x = 'Number of mutations',
      y = 'Cluster',
      subtitle = paste0("Mutational burden = ", sum(cl_tab_anno$nMuts), ', ', sum(cl_CCF$nDrivers), ' driver(s) in total.')
    ) +
    guides(
      fill = guide_legend('Number of drivers'),
      color = guide_legend(paste0('p < ', alpha_level))
    )
}

# Auxiliary function that makes the test used to color bars in plot_data_clone_size
pval_subcsz = function(x, p)
{
  cl = CCF_clusters(x, p)
  
  cltr = cl %>%
    filter(!is.clonal, !is.driver) %>%
    summarise(
      mu = mean(nMuts),
      sigma = sd(nMuts)
    )
  
  pv = function(y) { 1 - pnorm(y, mean = cltr$mu, sd = cltr$sigma) }
  
  cl %>%
    filter(!is.clonal, is.driver) %>%
    mutate(
      pvalue = pv(nMuts)
    ) %>%
    select(cluster, pvalue)
}


