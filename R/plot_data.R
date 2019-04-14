plot_data = function(x, p)
{
  sm = Samples(x, p)
  cl = CCF_clusters(x, p)
  st = Stats(x)

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
  cl_tab_anno$cluster = factor(cl_tab_anno$cluster, levels = cluster_ordering)
  cl_CCF$cluster = factor(cl_CCF$cluster, levels = cluster_ordering)

  # Drivers number
  mD = max(cl_CCF$nDrivers)
  cl_tab_anno$nDrivers = factor(cl_tab_anno$nDrivers, levels = paste(0:mD))

  pdata = ggplot(
    cl_tab,
    aes(
      x = region,
      y = cluster,
      z = CCF,
      fill = CCF,
      color = is.driver)
  ) +
    geom_tile(aes(width = .8, height = .8), size = 1) +
    geom_text(aes(label = CCF), color = 'orange') +
    scale_fill_distiller(palette = 'Blues', na.value = 'gainsboro', direction = 1) +
    scale_color_manual(
      values = c(`TRUE` = 'red', `FALSE` = NA)
    ) +
    theme_minimal(base_size = 10) +
    theme(
      legend.position = 'bottom',
      legend.key.size = unit(3, 'mm')
    ) +
    labs(
      title = paste("Clones for", p),
      subtitle = paste("Clonal cluster = ", cl$cluster[cl$is.clonal])
        ) +
    guides(fill = guide_colourbar(barwidth = 6))

  pvals = pval_subcsz(x, p)

  cl_tab_anno = cl_tab_anno %>%
    left_join(pvals, by = 'cluster') %>%
    mutate(
      # label = ifelse(is.na(pvalue), label, paste0(label, ' (p = ', round(pvalue, 3), ')')),
      significant = pvalue < 0.05
    )
  cl_tab_anno$cluster = factor(cl_tab_anno$cluster, levels = cluster_ordering)

  pmuts = ggplot(
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
    scale_color_manual(values = c(`TRUE` = 'forestgreen', `FALSE` = 'red')) +
    geom_text(aes(label = label), color = 'black', size = 2.5, hjust = 0) +
    # scale_y_discrete(sec.axis = dup_axis()) +
    theme_minimal(base_size = 10) +
    theme(
      legend.position = 'bottom',
      legend.key.size = unit(3, 'mm')
    ) +
    labs(
      title = "",
      subtitle = paste("Mutational burden = ", sum(cl_tab_anno$nMuts), ' - ', sum(cl_CCF$nDrivers), 'drivers')
    )

  # Clonal Subclonal
  pcoho = ggplot(
    st,
    aes(
      x = numTruncalMutations,
      y = numSubclonalMutations,
      color = numDriverMutations
    )
  ) +
    stat_density_2d(alpha = .3) +
    geom_point() +
    ggrepel::geom_label_repel(
      data = st %>% filter(patientID == !!p),
      aes(label = patientID),
      nudge_x = .2,
      nudge_y = .2
    ) +
    coord_cartesian(clip = 'off') +
    scale_color_gradient(low = 'steelblue', high = 'darkred') +
    theme_minimal(base_size = 10) +
    theme(
      legend.position = 'bottom',
      legend.key.size = unit(3, 'mm')
    ) +
    labs(
      title = paste("Cohort comparison"),
      subtitle = paste("Clonal versus subclonal mutations")
    ) +
    scale_x_log10() +
    scale_y_log10() +
    guides(color = guide_colourbar(barwidth = 6))

  ggpubr::ggarrange(
    pdata,
    pmuts,
    pcoho,
    ncol = 3,
    nrow = 1
  )

}

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
