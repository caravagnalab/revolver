require(tidyverse)

# Load your data
input = readr::read_tsv("~/Downloads/example_data.txt") %>% 
  mutate(cluster = paste(cluster))

require(revolver)

my_cohort = revolver_cohort(
  input %>% filter(!(patientID %in% c("S3", "S4"))), 
  CCF_parser = CCF_parser,
  ONLY.DRIVER = FALSE, 
  MIN.CLUSTER.SIZE = 5, # remove small clusters (some have just 1 mutation)
  annotation = "Test dataset"
)

# First, a print reports that  
# - Some driver variantIDs occur only once and should therefore be removed. 
print(my_cohort)

# This shows which one give the error
Stats_drivers(my_cohort) %>% 
  filter(N_tot == 1) 

# We remove those by using the variantID
my_cohort = remove_drivers(
  my_cohort,
  Stats_drivers(my_cohort) %>% 
    filter(N_tot == 1) %>% 
    pull(variantID)
)
  
# Attempt to compute trees
my_cohort = compute_clone_trees(my_cohort, sspace.cutoff = 1000, n.sampling = 500)

# The above command crashed at patient S7, raising errors
# [easypar] run 5 - Error in ClonEvol_surrogate(clusters, samples, clonal.cluster, min.CCF = 0.01): 
#   This patient has no trees, raising an error. Check you CCF estimates ...
# which means that S7 has CCF that are inconsistent with all the possible trees we can build.
# You should check your clustering.

# We remove that patient, this can remove further drivers
my_cohort = remove_patients(my_cohort, "S7")

# This gives no errors now, we can rebuild the trees
print(my_cohort)

my_cohort = compute_clone_trees(my_cohort, sspace.cutoff = 1000, n.sampling = 500)

# Trees are now there (see Trees per patient    : YES )
print(my_cohort)

# We can fit the cohort
my_cohort = revolver_fit(
  my_cohort, 
  parallel = F, 
  n = 3, 
  initial.solution = NA)

# .. compute clusters
my_cohort = revolver_cluster(
  my_cohort, 
  split.method = 'cutreeHybrid',
  min.group.size = 3)

plot_clusters(my_cohort, cutoff_trajectories = 1, cutoff_drivers = 0)
plot_drivers_graph(my_cohort)
plot_dendrogram(my_cohort)
plot_DET_index(my_cohort)
plot_drivers_clonality(my_cohort)
plot_patient_trees(my_cohort, "S5")

