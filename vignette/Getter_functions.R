data("TRACERx")

LF = list.files(path = '../R', all.files = TRUE, full.names = TRUE)
sapply(LF, source)

TRACERx = TRACERx %>% filter(patientID %in% c('CRUK0001', "CRUK0002", 'CRUK0003', 'CRUK0004', 'CRUK0005',
                                              'CRUK0006', 'CRUK0007', 'CRUK0008', 'CRUK0009',  'CRUK0010'))


# This will generate an erorr
# cohort = revolver_cohort(TRACERx)

# This will fix the input as required
TRACERx$cluster = paste(TRACERx$cluster)

# =-=-=-=-=-=-=-=-=-=-=-=-=-
# REVOLVER analysis
# =-=-=-=-=-=-=-=-=-=-=-=-=-

# Create the cohort
cohort = revolver_cohort(TRACERx,
                       options = list(ONLY.DRIVER = FALSE, MIN.CLUSTER.SIZE = 10))

# Create also phylogenetic trees
for(p in cohort$patients) cohort = revolver_compute_phylogenies(cohort, patient = p)

# Fit the REVOLVER model
cohort = revolver_fit(cohort, initial.solution = NA, n = 3)


# Simple dump of all the cohort
# all_cohort_plots = lapply(cohort$patients, plot_data, x = cohort)
# pdf("Cohort.pdf", width = 10, height = 5)
# lapply(all_cohort_plots, print)
# dev.off()

# =-=-=-=-=-=-=-=-=-=-=-=-=-
# Getter function for the tree
# =-=-=-=-=-=-=-=-=-=-=-=-=-
Phylo(cohort, 'CRUK0001')

# Top-rank tree
plot_trees(cohort, 'CRUK0001', rank = 1)
plot_trees(cohort, 'CRUK0001', rank = 1, icon = TRUE)

# Top-4 trees, arranged 2x2
ggarrange(
  plotlist = plot_trees(cohort, 'CRUK0001', rank = 1:4, icon = TRUE),
  nrow = 2, ncol = 2,
  labels = 1:4
  )

# =-=-=-=-=-=-=-=-=-=-=-=-=-
# Getter function for the fits
# =-=-=-=-=-=-=-=-=-=-=-=-=-
Fit(cohort, 'CRUK0001')

# =-=-=-=-=-=-=-=-=-=-=-=-=-
# Summary statistics for the cohort
# =-=-=-=-=-=-=-=-=-=-=-=-=-

# Cohort infos
Stats(cohort)
Stats(cohort, patients = "CRUK0001")

# Drivers
Stats_drivers(cohort)
Stats_drivers(cohort, drivers = c('APC', 'KRAS'))

# Trees
Stats_trees(cohort)
Stats_trees(cohort, patients = "CRUK0001")


# =-=-=-=-=-=-=-=-=-=-=-=-=-
# Fitting of a model with REVOLVER
# =-=-=-=-=-=-=-=-=-=-=-=-=-

cohort = revolver_fit(cohort,
                      initial.solution = NA)


# =-=-=-=-=-=-=-=-=-=-=-=-=-
# An example nice plot of data and tree for a patient
# =-=-=-=-=-=-=-=-=-=-=-=-=-

# Data plot
data_plot = plot_data(cohort, patient = 'CRUK0001')

# Trees pool icon format of 9 trees
icon_trees = ggarrange(
  plotlist = plot_trees(cohort, patient = 'CRUK0001', rank = 1:9, icon = TRUE),
  labels = 1:9,
  ncol = 3,
  nrow = 3)

# High-resolution top-rank plot
top_rank = plot_trees(cohort, 'CRUK0001', rank = 1)[[1]]

# Assembly
ggarrange(
  data_plot, # 1 x 3
  ggarrange(top_rank, icon_trees, # 1 x 2 scaled in width
            nrow = 1, ncol = 2,
            widths = c(1, 2),
            labels = c("Top rank", "")
  ),
  nrow = 2, # Overall figure stacked
  ncol = 1,
  labels = c('Data', '', '')
)
