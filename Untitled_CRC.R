load('../revolver copy master 12 otto/data_old/CRC.RData')
CRC

require(tidyverse)
require(revolver)

CROSS_CRC_ADENOCARCINOMA = CRC %>% as_tibble %>% mutate(cluster = paste(cluster))

save(CROSS_CRC_ADENOCARCINOMA, file = "~/Downloads/CROSS_CRC_ADENOCARCINOMA.RData")

CROSS_CRC_ADENOCARCINOMA_REVOLVER = revolver_cohort(CROSS_CRC_ADENOCARCINOMA, MIN.CLUSTER.SIZE = 0)

non_recurrent = Stats_drivers(CROSS_CRC_ADENOCARCINOMA_REVOLVER) %>% 
  filter(N_tot == 1) %>% 
  pull(variantID)

CROSS_CRC_ADENOCARCINOMA_REVOLVER = remove_drivers(CROSS_CRC_ADENOCARCINOMA_REVOLVER, non_recurrent)

CROSS_CRC_ADENOCARCINOMA_REVOLVER = compute_mutation_trees(CROSS_CRC_ADENOCARCINOMA_REVOLVER)

CROSS_CRC_ADENOCARCINOMA_REVOLVER = revolver_fit(CROSS_CRC_ADENOCARCINOMA_REVOLVER, parallel = F, n = 5, initial.solution = NA)

CROSS_CRC_ADENOCARCINOMA_REVOLVER = revolver_cluster(CROSS_CRC_ADENOCARCINOMA_REVOLVER, 
                                                     split.method = 'cutreeHybrid',
                                                     min.group.size = 3)

CROSS_CRC_ADENOCARCINOMA_REVOLVER = revolver_jackknife(CROSS_CRC_ADENOCARCINOMA_REVOLVER, 
                                                       resamples = 10, 
                                                       ... = list(
                                                         `split.method` = 'cutreeHybrid',
                                                         `min.group.size` = 3
                                                         ) 
                                                       )

plot_clusters(CROSS_CRC_ADENOCARCINOMA_REVOLVER, cutoff_trajectories = 1, cutoff_drivers = 0)

plot_drivers_graph(CROSS_CRC_ADENOCARCINOMA_REVOLVER)
plot_dendrogram(CROSS_CRC_ADENOCARCINOMA_REVOLVER)
plot_DET_index(CROSS_CRC_ADENOCARCINOMA_REVOLVER)
plot(CROSS_CRC_ADENOCARCINOMA_REVOLVER)
plot_drivers_clonality(CROSS_CRC_ADENOCARCINOMA_REVOLVER)
plot_drivers_occurrence(CROSS_CRC_ADENOCARCINOMA_REVOLVER)
plot_jackknife_cluster_stability(CROSS_CRC_ADENOCARCINOMA_REVOLVER)
plot_jackknife_coclustering(CROSS_CRC_ADENOCARCINOMA_REVOLVER)
plot_jackknife_trajectories_stability(CROSS_CRC_ADENOCARCINOMA_REVOLVER)
plot_patient_trees(CROSS_CRC_ADENOCARCINOMA_REVOLVER, CROSS_CRC_ADENOCARCINOMA_REVOLVER$patients[1])


# Breast fits

load('../revolver copy master 12 otto/data_old/Breast.fit.RData')
YATES_BREAST_CANCERS = Breast.fit$dataset %>% as_tibble

YATES_BREAST_CANCERS = YATES_BREAST_CANCERS %>% as_tibble %>% mutate(cluster = paste(cluster))

save(YATES_BREAST_CANCERS, file = "~/Downloads/YATES_BREAST_CANCERS.RData")

YATES_BREAST_CANCERS_REVOLVER = revolver_cohort(YATES_BREAST_CANCERS, MIN.CLUSTER.SIZE = 0)

YATES_BREAST_CANCERS_REVOLVER = compute_mutation_trees(YATES_BREAST_CANCERS_REVOLVER)

YATES_BREAST_CANCERS_REVOLVER = revolver_fit(YATES_BREAST_CANCERS_REVOLVER, parallel = F, n = 5, initial.solution = NA)

YATES_BREAST_CANCERS_REVOLVER = revolver_cluster(YATES_BREAST_CANCERS_REVOLVER, 
                                                     split.method = 'cutreeHybrid',
                                                     min.group.size = 3)

YATES_BREAST_CANCERS_REVOLVER = revolver_jackknife(YATES_BREAST_CANCERS_REVOLVER, 
                                                       resamples = 10, 
                                                       ... = list(
                                                         `split.method` = 'cutreeHybrid',
                                                         `min.group.size` = 3
                                                       ) 
)

plot_clusters(YATES_BREAST_CANCERS_REVOLVER, cutoff_trajectories = 1, cutoff_drivers = 0)

plot_drivers_graph(YATES_BREAST_CANCERS_REVOLVER)
plot_dendrogram(YATES_BREAST_CANCERS_REVOLVER)
plot_DET_index(YATES_BREAST_CANCERS_REVOLVER)
plot(YATES_BREAST_CANCERS_REVOLVER)
plot_drivers_clonality(YATES_BREAST_CANCERS_REVOLVER)
plot_drivers_occurrence(YATES_BREAST_CANCERS_REVOLVER)
plot_jackknife_cluster_stability(YATES_BREAST_CANCERS_REVOLVER)
plot_jackknife_coclustering(YATES_BREAST_CANCERS_REVOLVER)
plot_jackknife_trajectories_stability(YATES_BREAST_CANCERS_REVOLVER)
plot_patient_trees(YATES_BREAST_CANCERS_REVOLVER, YATES_BREAST_CANCERS_REVOLVER$patients[1])
plot_patient_oncoprint(YATES_BREAST_CANCERS_REVOLVER, YATES_BREAST_CANCERS_REVOLVER$patients[1])
