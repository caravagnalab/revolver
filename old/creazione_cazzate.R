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

TRACERx_cohort = revolver_compute_phylogenies(TRACERx_cohort, "CRUK0002", 
                                              options = list(sspace.cutoff = 10000,
                                                                n.sampling = 5000,
                                                                overwrite = TRUE,
                                                                store.max = 100)
                                              )
TRACERx_cohort = revolver_compute_phylogenies(TRACERx_cohort, "CRUK0001")

load('data_old/CRC.RData')

CRC_cohort = revolver_cohort(
  CRC %>% as_tibble() %>% mutate(cluster = paste(cluster)),
  CCF_parser = revolver::CCF_parser,
  annotation =   "Cross et al, colorectal adenocarcinoma.", 
  options = list(ONLY.DRIVER = FALSE, MIN.CLUSTER.SIZE = 0)
)

CRC_cohort = revolver_compute_CLtrees(CRC_cohort, CRC_cohort$patients[3], 
                                              options = list(sspace.cutoff = 10000,
                                                             n.sampling = 5000,
                                                             overwrite = TRUE,
                                                             store.max = 100)
)
