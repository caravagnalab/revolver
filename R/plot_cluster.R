#' Plot the heatmaps of REVOLVER"s clusters.
#'
#' @description
#'
#' Plot the heatmaps of REVOLVER"s clusters, as tiles.
#'
#' The top tile is patients vs trajectories, and bottom is patients vs drivers.
#' For drivers colours reflect mean CCF/ binary values of a driver in every patient, and
#' clonality status. For trajectories colours reflect if they are initiating or progressing,
#' depending on the present of GL in the trajectory.
#'
#' Patients are sorted by cluster to match the dendrogram that one can obtain with \code{\link{plot_dendrogram}}.
#'
#' @param x A \code{REVOLVER} object with fits and clusters.
#' @param clusters_palette A palette function that should return the colour of
#' an arbitrary number of clusters.
#' @param cutoff_drivers Plot only drivers that occur in at least \code{cutoff_drivers} patients.
#' @param cutoff_trajectories Plot only trajectories that occur in at least \code{cutoff_trajectories} patients.
#' @param arrow.symbol UNICODE code to display arrows. Saving to PDF outputs with standard
#' methods (ggsave, cairo, etc), often can lead to errors with UNICODE chars; therefore either print to PNG or change
#' this variable to, e.g.,  \code{" --> "} to render arrows properly.
#'
#' @family Plotting functions
#'
#' @return A \code{ggplot} plot.
#' @export
#'
#' @examples
#' # Data released in the 'evoverse.datasets'
#' data('TRACERx_NEJM_2017_REVOLVER', package = 'evoverse.datasets')
#'
#' plot_clusters(TRACERx_NEJM_2017_REVOLVER)
#'
#' plot_clusters(TRACERx_NEJM_2017_REVOLVER, cutoff_drivers = 10, cutoff_trajectories = 3)
plot_clusters = function(x,
                         cluster_palette = distinct_palette_few,
                         cutoff_drivers = 5,
                         cutoff_trajectories = 4,
                         arrow.symbol = ' \u2192 ')
{
  obj_has_clusters(x)

  # features = get_features(x, patients)
  features = get_features(x)

  # =-=-=-=-=-=-=-=-
  # Prepare driver and trajectory data required for this plot
  # =-=-=-=-=-=-=-=-

  # Drivers occurring in at least cutoff_drivers patients
  occurrence_drivers = reshape2::melt(features$Matrix_drivers, id = 'patientID') %>%
    as_tibble %>%
    dplyr::rename(driverID = variable, n = value) %>%
    dplyr::group_by(driverID) %>%
    dplyr::summarise(n = sum(n)) %>%
    dplyr::arrange(desc(n)) %>%
    dplyr::filter(n > cutoff_drivers) %>%
    dplyr::pull(driverID) %>%
    as.character

  # Mean CCF for occurrence_drivers, augmented with clonality status for the driver event
  mean_CCF = reshape2::melt(features$Matrix_mean_CCF, id = 'patientID') %>%
    as_tibble %>%
    dplyr::rename(driverID = variable, mean_CCF = value) %>%
    dplyr::filter(driverID %in% occurrence_drivers) %>%
    dplyr::filter(mean_CCF > 0)

  clonal_events = reshape2::melt(features$Matrix_clonal_drivers, id = 'patientID') %>%
    as_tibble %>%
    dplyr::rename(driverID = variable, clonal = value)

  mean_CCF = mean_CCF %>%
    dplyr::full_join(clonal_events, by = c('patientID', 'driverID')) %>%
    dplyr::mutate(
      clonal = ifelse(clonal == 1, "clonal", 'subclonal'),
      driverID = paste(driverID)
    )

  # Trajectories in at least cutoff_trajectories
  occurrence_trajectories = reshape2::melt(features$Matrix_trajectories, id = 'patientID') %>%
    as_tibble %>%
    dplyr::rename(trajectory = variable, n = value) %>%
    dplyr::group_by(trajectory) %>%
    dplyr::summarise(n = sum(n)) %>%
    dplyr::arrange(desc(n)) %>%
    dplyr::filter(n > cutoff_trajectories) %>%
    dplyr::pull(trajectory) %>%
    as.character

  # Map of the required trajectories as of occurrence_trajectories
  map_occurrence_trajectories = reshape2::melt(features$Matrix_trajectories, id = 'patientID') %>%
    as_tibble %>%
    dplyr::rename(trajectory = variable, n = value) %>%
    dplyr::filter(trajectory %in% occurrence_trajectories) %>%
    dplyr::filter(n > 0) %>%
    tidyr::separate(
      col = trajectory,
      into = c('from', 'to'),
      sep = ' --> ',
      remove = FALSE
    ) %>%
    dplyr::mutate(trajectory = paste(trajectory))

  # Add also clusters to trajectories
  map_occurrence_trajectories = map_occurrence_trajectories %>%
    dplyr::mutate(trajectory = paste(trajectory)) %>%
    dplyr::left_join(Cluster(x),
              by = 'patientID')

  # =-=-=-=-=-=-=-=-
  # Work out orderings that will be used to make everything consistent across plots,
  # special unicode symbols and colors for nicer plots
  # =-=-=-=-=-=-=-=-

  # Dendrogram - it gives the ordering of the patients which are displayed on the x-axcis
  hc = x$cluster$fits$hc

  # Patient ordering - from the dedrogram, thi defines the levels of the factors used
  factors_patient_level = hc$order.lab
  
  cl_ordering = Cluster(x) %>% 
    arrange(cluster) %>% 
    pull(patientID)
  
  # if(any(cl_ordering != factors_patient_level))
  # {
  #   warning("Dendrogram ordering does not reflect clustering - strange behaviour of the cutting algorithms?")
  #   # factors_patient_level = cl_ordering
  # }

  # Get colors for the clusters
  clusters_colors = get_cluster_colors(x, cluster_palette)

  # Assign the colors following factors_patient_level
  patients_factors_colors = sapply(factors_patient_level,
                                   function(y)
                                     clusters_colors[Cluster(x, y) %>% dplyr::pull(cluster)])
  names(patients_factors_colors) = factors_patient_level

  bars_separation = Cluster(x, factors_patient_level)
  bars_separation$cluster = factor(bars_separation$cluster, levels = unique(bars_separation$cluster))
  bars_separation = bars_separation %>% pull(cluster) %>% table %>% cumsum + 0.5

  # Trajectories ordering - by overall frequency
  factors_trajectory_level = map_occurrence_trajectories %>%
    dplyr::group_by(trajectory) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::arrange(desc(n)) %>%
    dplyr::pull(trajectory) %>%
    rev

  # Drivers ordering - by overall frequency
  factors_drivers_level = occurrence_drivers %>% rev

  # Add unicode char for -->
  factors_trajectory_level = gsub(factors_trajectory_level,
                                  pattern = ' --> ',
                                  replacement = arrow.symbol)
  map_occurrence_trajectories$trajectory = gsub(
    map_occurrence_trajectories$trajectory,
    pattern = ' --> ',
    replacement = arrow.symbol
  )

  # Number of clusters
  nclusters = Cluster(x) %>% pull(cluster) %>% unique %>% length

  # Ggplot trajectories
  trajectories_plot = ggplot(
    map_occurrence_trajectories %>%
      mutate(type = ifelse(from == 'GL', "Initiation", "Progression")),
    aes(x = patientID, y = trajectory)
  ) +
    geom_tile(aes(
      width = .8,
      height = .8,
      fill = type
    ), size = .5) +
    geom_vline(
      xintercept = bars_separation,
      size = .3,
      color = 'darkred',
      linetype = 'dashed'
    ) +
    my_ggplot_theme() +
    theme(axis.text.x = element_blank(),
          legend.position = 'top') +
    labs(
      x = NULL,
      y = 'Trajectory',
      title = 'REVOLVER trajectories and drivers',
      subtitle = x$annotation
    ) +
    scale_fill_manual(values = c(`Initiation` = 'darkred', `Progression` = 'steelblue')) +
    guides(fill = guide_legend("Trajectory type")) +
    scale_color_manual(values = clusters_colors) +
    scale_x_discrete(limits = factors_patient_level) +
    scale_y_discrete(limits = factors_trajectory_level)

  # Ggplot drivers
  # pt_ord = mean_CCF %>% filter(mean_CCF > 0) %>% pull(patientID)
  
  driver_plot =
    ggplot(
      mean_CCF %>% filter(mean_CCF > 0),
      aes(
        x = patientID,
        y = driverID,
        fill = mean_CCF,
        color = clonal
      )
    ) +
    geom_tile(aes(width = .8, height = .8), size = .5) +
    labs(x = 'Patient',
         y = 'Driver',
         caption = paste0('k = ', nclusters, ' clusters, ',
                          'n = ', x$n$patients, ' patients.')
         ) +
    my_ggplot_theme() +
    theme(axis.text.x = element_text(
      angle = 90,
      size = 5
      # color = patients_factors_colors[pt_ord]
    )) +
    guides(
      fill = guide_colorbar("Mean across biopsies  ", barwidth = unit(3, 'cm')),
      color = guide_legend("Clonality status")
    ) +
    # scale_fill_viridis_c(option = 'D', direction = -1, values = seq(0, 1, 0.1), limits = c(0, 1)) +
    scale_fill_distiller(palette = 'Blues',
                         direction = 1,
                         limits = c(0, 1)) +
    scale_color_manual(values = c(`clonal` = 'darkgoldenrod2', `subclonal` = NA)) +
    scale_x_discrete(limits = factors_patient_level) +
    scale_y_discrete(limits = factors_drivers_level) +
    geom_vline(
      xintercept = bars_separation,
      size = .3,
      color = 'darkred',
      linetype = 'dashed'
    )

  clbar = Cluster(x) %>%
    ggplot(aes(x = patientID, y = 1)) +
    geom_tile(aes(fill = cluster)) +
    my_ggplot_theme() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x  = element_blank(),
      axis.title.x = element_blank()
      # axis.text.x = element_text(angle = 90, size = 5)
    ) +
    labs(y = '', x = NULL) +
    scale_fill_manual(values = clusters_colors) +
    guides(fill = 'none') +
    scale_x_discrete(limits = factors_patient_level) 
    
  cowplot::plot_grid(
    trajectories_plot,
    clbar,
    driver_plot,
    nrow = 3,
    ncol = 1,
    rel_heights = c(1, .1, 1),
    align = 'v'
  )
}
