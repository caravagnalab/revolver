load('data/TRACERx_data.rda')
subset_data = TRACERx_data %>% 
  filter(patientID %in% paste0('CRUK000', 1:9))

x = revolver_cohort(subset_data, ONLY.DRIVER = T, MIN.CLUSTER.SIZE = 0)

Stats_cohort(x)
Stats_drivers(x)
Stats_trees(x)

library(ctree)

x = compute_clone_trees(x, patients = x$patients, sspace.cutoff=1000, n.sampling = 200, store.max = 100)

data(TRACER_NEJM_2017, package = 'REVOLVER.datasets')

ctree::plot_CCF_clusters(Phylo(x, "CRUK0002", rank = 1))
plot_CCF_histogram(x, 'CRUK0002')
plot_drivers_clonality(x)
plot_drivers_occurrence(x)
plot_patient_data(x, 'CRUK0002')

x = revolver_fit(x)
x = revolver_cluster(x)
plot_clusters(x)
 
Stats_fits(x)
