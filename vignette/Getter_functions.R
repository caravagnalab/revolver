data("TRACERx")

# This will generate an erorr
cohort = revolver_cohort(TRACERx)

# This will fix the input as required
TRACERx$cluster = paste(TRACERx$cluster)

cohort = revolver_cohort(TRACERx,
                       options = list(ONLY.DRIVER = FALSE, MIN.CLUSTER.SIZE = 10))

Data(cohort, 'CRUK0001')
Drivers(cohort, 'CRUK0001')

Samples(cohort, 'CRUK0001')

Truncal(cohort, 'CRUK0001')
Subclonal(cohort, 'CRUK0001')

CCF(cohort, 'CRUK0001')
CCF_clusters(cohort, 'CRUK0001')

Phylo(cohort, 'CRUK0001')
Fit(cohort, 'CRUK0001')

Stats(cohort)

plot_data(cohort, 'CRUK0001')
plot_data(cohort, 'CRUK0062')
plot_data(cohort, 'CRUK0069')

all_cohort_plots = lapply(cohort$patients, plot_data, x = cohort)
pdf("Cohort.pdf", width = 10, height = 5)
lapply(all_cohort_plots, print)
dev.off()


