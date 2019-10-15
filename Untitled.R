data('TRACERx_NEJM_2017_REVOLVER', package = 'evoverse.datasets')

x = remove_patients(x = TRACERx_NEJM_2017_REVOLVER, patientID = paste0("CRUK000", 1:9))


x = revolver::revolver_jackknife(TRACERx_NEJM_2017_REVOLVER, resamples = 2, parallel = F)
