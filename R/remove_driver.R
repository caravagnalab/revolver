
#' Remove a driver event from the cohort.
#'
#' @description Each event is identied through its ID (variantID); 
#' with this function, you can remove a driver events, which consists
#' in flagging it as FALSE in the \code{is.Driver} column of the data, and updating
#' the information transfer. If you have fit the models, however, you should 
#' re-run the fit after this modification because this modification does not propagate.
#'
#' @param x A REVOLVER cohort.
#' @param variantID ID of the driver event to remove.
#'
#' @return A modified cohort where the required event is no longer a driver.
#' @export
#'
#' @examples
#' data(CRC.cohort)
#' head(CRC.cohort$dataset) # Info to use to remove
#' new.data = revolver_removeDriver(CRC.cohort, "adenoma_1", "APC", "NOTHING", "1")
#' head(new.data$dataset)
remove_driver = function(
  x,
  variantID
  )
{
  pioStr('Removing driver event', variantID, suffix = '\n')
  Stats_drivers(x, drivers = variantID) %>% pioDisp
  
  N = length(x$patients)
  
  pb =  dplyr::progress_estimated(n = N, min_time = 2)
  
  for(patient in x$patients)
  {
    pb$tick()$print()
    
    # Check if the patient has this driver, skip if not
    drv = Drivers(x, patient) %>% pull(variantID)
    
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
  
  # if has fits, suggest to recompute
  if(has_fits(x)) {
    message("The cohort has fits, changing drivers require re-computing the fits. 
            Please re-run revolver_fit")
  }
  
  if(has_clusters(x)) {
    message("The cohort has clusters, changing drivers require re-computing the clustering. 
            Please re-run revolver_cluster")
  }
  
  return(x)
}