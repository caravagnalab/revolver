#' Remove a list of driver events from the cohort.
#'
#' @description Each event is identified through its id (\code{variantID});
#' with this function, you can remove a list of driver events, which consists
#' in flagging them as \code{FALSE} in the \code{is.Driver} column of the data, and updating
#' the information transfer. If you have fit the models or clustered the cohort,
#' you should re-run the analyses after this modification; for this reason,
#' any previous result from those analyses is cancelled from the returned object.
#'
#' Notice also that some patients might be removed by this function, because if
#' they have no longer drivers then they cannot be fit afterwards.
#'
#' @param x A REVOLVER cohort.
#' @param variantID Id of the driver event to remove.
#'
#' @return A modified cohort where the required events are no longer annotated as drivers.
#' @export
#'
#' @examples
#' # Data released in the 'evoverse.datasets'
#' data('TRACERx_NEJM_2017_REVOLVER', package = 'evoverse.datasets')
#'
#' print(TRACERx_NEJM_2017_REVOLVER)
#'
#' new_cohort = remove_drivers(TRACERx_NEJM_2017_REVOLVER, "MET")
#' print(new_cohort)
remove_drivers = function(
  x,
  variantID,
  check = TRUE
  )
{
  if(!all(variantID %in% x$variantIDs.driver))
    stop(variantID, "is not a driver in the cohort, aborting.")

  cli::cli_rule('Removing driver events', right = paste(variantID, collapse = ', '))
  cat("\n")

  Stats_drivers(x, drivers = variantID) %>% pioDisp

  N = length(x$patients)

  pb =  dplyr::progress_estimated(n = N, min_time = 2)

  # List of patients that stay after this editing
  new_patients = c()

  for(patient in x$patients)
  {
    pb$tick()$print()

    # Drivers for this patient
    drv = Drivers(x, patient) %>% pull(variantID)

    # Check  if the patient has other drivers, because if it has not then the
    # patient needs to be removed. There cannot be a patient without drivers in REVOLVER
    # (also the tree-based packages would throw an error)
    if(all(drv %in% variantID)) 
    {
      warning(patient, " will have no drivers, and therefore will be removed.\n")

      next;
    }

    # Then at this point the patient stays in, so we add it to list that tracks patients
    new_patients = c(new_patients, patient)

    # If the patient has not this driver, we skip it
    if(!any(variantID %in% drv)) next;

    # Remove from the data of the patient
    new_data = Data(x, patient) %>%
      mutate(is.driver = ifelse(variantID %in% !!variantID, FALSE, is.driver))

    CCFclone_headers = new_data %>%
      group_by(cluster) %>%
      summarise(
        nMuts = n(),
        is.driver = any(is.driver),
        is.clonal = all(is.clonal)
      ) %>%
      arrange(desc(nMuts))
    
    CCFclone_values = new_data %>%
      select(cluster, Samples(x, patient)) %>%
      reshape2::melt(id = 'cluster') %>%
      group_by(cluster, variable) %>%
      summarise(CCF = median(as.numeric(value))) %>%
      spread(variable, CCF) %>%
      ungroup()
    
    new_CCF_clusters = CCFclone_headers %>% full_join(CCFclone_values, by = 'cluster')

    # Clones with still a driver
    # clones_driver_profile = new_data %>%
    #   group_by(cluster) %>%
    #   summarise(is.driver = any(is.driver)) %>%
    #   filter(is.driver) %>%
    #   pull(cluster)
    # 
    # # Where the driver maps to
    # clones_assignment = new_data %>%
    #   filter(variantID %in% !!variantID) %>%
    #   pull(cluster)
    # 
    # # In this case, removing this driver we have no longer
    # # a driver clone, need to update the CCF data
    # if(!(clones_assignment %in% clones_driver_profile))
    #   new_CCF_clusters = new_CCF_clusters %>%
    #     mutate(is.driver = ifelse(cluster == clones_assignment, FALSE, is.driver))

    # Store the new data and CCF clusters
    x$dataset[[patient]] = new_data
    x$CCF[[patient]] = new_CCF_clusters

    # If it has the trees, we need to update them
    if(has_patient_trees(x, patient))
    {
      trees = Phylo(x, patient)

      # Get matrices from these objects, and remove the GL columns/ rows
      matrices = lapply(trees, function(x) x$adj_mat[rownames(x$adj_mat) != 'GL', colnames(x$adj_mat) != 'GL'])

      # Get scores for these matrices - vector format
      scores = sapply(trees, function(x) x$score)

      # Add these trees x$phylogenies[[patient]] =
      new_trees = easypar::run(
        FUN = function(w){
          ctree::ctree(
            new_CCF_clusters,
            new_data %>% filter(is.driver),
            Samples(x, patient),
            patient = patient,
            M = matrices[[w]],
            score = scores[[w]]
          )
        },
        PARAMS = lapply(seq_along(trees), list),
        parallel = FALSE,
        progress_bar = FALSE
      )

      if(easypar::numErrors(new_trees) > 0)
        stop("Errors updating the trees of ", patient, ' ~ Aborting.')

      x$phylogenies[[patient]] = new_trees
    }
  }

  # Retianed patients
  cli::cli_alert_info("Retained {.field {length(new_patients)}} patients after driver removal..\n")

  # if has fits, force to recompute
  if(has_fits(x)) {
    message("- The cohort has fits which will be cancelled, please re-compute the fits ...")
    x$fit = NULL
  }

  # if has fits, force to recompute
  if(has_clusters(x)) {
    message("- The cohort has clusters which will be cancelled, please re-compute the clusters ...")
    x$cluster = NULL
  }

  # Update counts and remove non-used patients
  x$patients = new_patients

  x$n$patients = length(x$patients)
  x$n$variants = x$n$variants - 1
  x$variantIDs.driver = setdiff(x$variantIDs.driver, variantID)
  x$n$drivers = x$n$drivers - 1

  x$dataset = x$dataset[new_patients]
  x$CCF = x$CCF[new_patients]
  x$phylogenies = x$phylogenies[new_patients]

  # Check the cohort
  if(check) revolver_check_cohort(x)

  return(x)
}
