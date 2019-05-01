TRACERx$cluster = paste(TRACERx$cluster)
TRACERx_data = TRACERx_data %>% as_tibble()

data(TRACERx_data)

usethis::use_data(TRACERx_data, overwrite = T)

TRACERx_cohort = revolver_cohort(
  TRACERx_data,
  CCF_parser = revolver::CCF_parser,
  annotation =   "TRACERx cohort released with the REVOLVER package"
)

usethis::use_data(TRACERx_cohort, overwrite = TRUE)

TRACERx_cohort = revolver_compute_phylogenies(TRACERx_cohort, "CRUK0002")
