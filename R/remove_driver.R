#' Remove a driver event from the cohort.
#'
#' @description Each event is identied through its id (\code{variantID});
#' with this function, you can remove a driver events, which consists
#' in flagging it as \code{FALSE} in the \code{is.Driver} column of the data, and updating
#' the information transfer. If you have fit the models or clustered the cohort,
#' you should re-run the analyses after this modification; for this reason,
#' any previous result from those analyses is cancelled from the returned object.
#'
#' Notice also that some patients might be removed by this function, because if
#' they have only one driver then they cannot be fit after driver removal.
#'
#' @param x A REVOLVER cohort.
#' @param variantID Id of the driver event to remove.
#'
#' @return A modified cohort where the required event is no longer a driver.
#' @export
#'
#' @examples
#' # Data released in the 'evoverse.datasets'
#' data('TRACERx_NEJM_2017_REVOLVER', package = 'evoverse.datasets')
#'
#' print(TRACERx_NEJM_2017_REVOLVER)
#'
#' new_cohort = remove_driver(TRACERx_NEJM_2017_REVOLVER, "MET")
#' print(new_cohort)
remove_driver = function(
  x,
  variantID
  )
{
  if(!(variantID %in% x$variantIDs.driver))
    stop(variantID, "is not a driver in the cohort, aborting.")

  pioStr('Removing driver event', variantID, suffix = '\n')
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
    if(length(drv) == 1 && drv == variantID) {
      warning(patient, " has only one driver (", variantID, ") and therefore will be removed.\n")

      next;
    }

    # Then at this point the patient stays in, so we add it to list that tracks patients
    new_patients = c(new_patients, patient)

    # If the patient has not this driver, we skip it
    if(!(variantID %in% drv)) next;

    # Remove from the data of the patient
    new_data = Data(x, patient) %>%
      mutate(is.driver = ifelse(variantID == !!variantID, FALSE, is.driver))

    new_CCF_clusters = CCF_clusters(x, patient)

    # Clones with still a driver
    clones_driver_profile = new_data %>%
      group_by(cluster) %>%
      summarise(is.driver = any(is.driver)) %>%
      filter(is.driver) %>%
      pull(cluster)

    # Where the driver maps to
    clones_assignment = new_data %>%
      filter(variantID == !!variantID) %>%
      pull(cluster)

    # In this case, removing this driver we have no longer
    # a driver clone, need to update the CCF data
    if(!(clones_assignment %in% clones_driver_profile))
      new_CCF_clusters = new_CCF_clusters %>%
        mutate(is.driver = ifelse(cluster == clones_assignment, FALSE, is.driver))

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
  pioStr("Retained patients", paste(new_patients, collapse = ', '), suffix = '\n')

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
  revolver_check_cohort(x)

  return(x)
}
