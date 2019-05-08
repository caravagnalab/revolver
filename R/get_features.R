#' Return summary features for the cohort.
#'
#' @description 
#' Computes a set of matrix-like summary features for the cohort, in the
#' form of a matrix. It will return the matrix of driver events with mean
#' CCF, their occurrence across all patients, and at the clonal and subclonal
#' level, and the occurrence of all evolutionary trajectories across patients
#' 
#' @param x 
#' @param patients 
#'
#' @return A set of matrices that report summary features of the models.
#' @export
#'
#' @examples
#' TODO
get_features = function(x, patients = x$patients)
{
  Np = length(patients)
  Nd = x$n$drivers
  
  # =-=-=-=-=-=-=-=-
  # Get all data that we need for driver calls across patients
  # =-=-=-=-=-=-=-=-
  All_drivers = lapply(patients, function(p) {
    samples = Samples(x, p)
    
    drivers = Drivers(x, p) 
    drivers$mean_CCF = apply(drivers[, samples], 1, mean)

    drivers
    })
  
  All_drivers = Reduce(bind_rows, All_drivers) %>%
    select(variantID, patientID, is.clonal, mean_CCF) %>%
    mutate(value = 1)
  
  # =-=-=-=-=-=-=-=-
  # Matrix of mean CCF
  # =-=-=-=-=-=-=-=-
  Matrix_mean_CCF = All_drivers %>% 
    select(variantID, patientID, mean_CCF) %>%
    spread(variantID, mean_CCF) %>%
    replace(is.na(.), 0)
  
  # =-=-=-=-=-=-=-=-
  # 0/1 drivers matrices, all mutations and then split by clonality status
  # =-=-=-=-=-=-=-=-
  
  # All together is trivial
  Matrix_drivers = All_drivers %>% 
    select(variantID, patientID, value) %>%
    spread(variantID, value) %>%
    replace(is.na(.), 0)
  
  # Only clonal ones
  Matrix_clonal_drivers = All_drivers %>% 
    filter(is.clonal) %>%
    select(variantID, patientID, value) %>%
    spread(variantID, value) %>%
    replace(is.na(.), 0)
  
  # Only subclonal
  Matrix_subclonal_drivers = All_drivers %>% 
    filter(!is.clonal) %>%
    select(variantID, patientID, value) %>%
    spread(variantID, value) %>%
    replace(is.na(.), 0)
  
  # Now we want to make them all the same dimension
  # to standardize this output... we create a template matrix and use it
  # to complement missing columns in each of the clonal/ subclonal ones
  complement = function(M, N) 
  {
    # missing patients and driver genes
    miss_Pat = setdiff(M$patientID, N$patientID)
    miss_drv = setdiff(colnames(M), colnames(N))
    
    # Add template 0-ed matrices with the right rows/ columns
    if(length(miss_Pat) > 0)
    {
      empty = M %>% filter(patientID %in% !!miss_Pat)
      empty[, 2:ncol(empty)] = 0
      
      N = bind_rows(N, empty)
    }
    
    if(length(miss_drv) > 0)
      N = bind_cols(N,
                    M %>% 
                      select(!!miss_drv) %>%
                      replace(TRUE, 0)
      )
     
    N[, colnames(M)]             
  }
  
  Matrix_clonal_drivers = complement(Matrix_drivers, Matrix_clonal_drivers) %>%
    replace(is.na(.), 0)
  Matrix_subclonal_drivers = complement(Matrix_drivers, Matrix_subclonal_drivers) %>%
    replace(is.na(.), 0)
  
  # =-=-=-=-=-=-=-=-
  # Get all data that we need for trajectories
  # =-=-=-=-=-=-=-=-
  All_trajectories = lapply(patients, function(p) {
    ITransfer(x, p, rank = 1, type = 'drivers', data = 'fits') %>%
      mutate(patientID = p)
  })
  
  All_trajectories = Reduce(bind_rows, All_trajectories)  %>%
    mutate(
      trajectory = paste0(from, ' --> ', to),
      value = 1
      ) %>%
    select(trajectory, patientID, value)

  Matrix_trajectories = All_trajectories %>% 
    spread(trajectory, value) %>%
    replace(is.na(.), 0)
  
  return(
    list(
      Matrix_mean_CCF = Matrix_mean_CCF,
      Matrix_drivers = Matrix_drivers,
      Matrix_clonal_drivers = Matrix_clonal_drivers,
      Matrix_subclonal_drivers = Matrix_subclonal_drivers,
      Matrix_trajectories = Matrix_trajectories
    )
  )
}
