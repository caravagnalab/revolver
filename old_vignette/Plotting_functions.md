Plotting functions
================
Giulio Caravagna
May 1, 2019

``` r
# Load REVOLVER 
library(revolver)
```

    ##  [ REVOLVER - Repeated Evolution in Cancer ] 
    ##  Author :  Giulio Caravagna <gcaravagn@gmail.com> 
    ##  GitHub :  caravagn/revolver

**Note** Since the new release of REVOLVER, we have implemented the
almost all plotting functions within the `ggplot` domain, using packages
like `ggraph`, `ggpubr`, `ggrepel` etc. This means that output plots can
often be assembled via standard methods (`ggpubr`, `cowplot`, etc.).

In this vignette we make use of one of the cohort objects released with
the tool.

``` r
# Data from the TRACERx cohort 
load("Cohort.RData") 

# Load a plot-assembly package
require(ggpubr)
```

# S3 plot objects for a cohort

S3 cohort objects before (`rev_cohort`) and after (`rev_cohort_fit`) the
fit have their own `plot` methods. For a cohort, the plot shows a
scatter of the patients’ mutatinal burden (clonal versus subclonal), the
drivers annotated and the number of clones annotated; this visualization
can be computed even if trees for the cohort have not yet been created.
After fit, this scatterplot changes showing the overall mutaitonal
burden versus one minus the penalty of the fit, which gives idea of
which patients have more homogenous evolutionary trajectories.

``` r
# The class of this object, 
class(cohort)
```

    ## [1] "rev_cohort_fit"

``` r
# We force call of both of them explicitely because rev_cohort_fit inherits from rev_cohort
ggarrange(
  revolver:::plot.rev_cohort(cohort),
  revolver:::plot.rev_cohort_fit(cohort),
  nrow = 1,
  ncol = 2
)
```

![](Plotting_functions_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

# Plotting data

You can plot the CCF values for each clone annotated in a patient’s
data, as well as the clone size for each clone annotated in a patient,
where subclones with drivers are annotated if they are significantly
larger than subclones without drivers. See `?plot_data_clone_size` for
information about how the significance is
computed.

``` r
# Here the subclone with NF1 and ARHGAP35 are significantly larger than the subclones without drivers.
# The subclone with PASK however is not.
ggarrange(
  plot_data_clusters(cohort, 'CRUK0001'),
  plot_data_clone_size(cohort, 'CRUK0001'),
  nrow = 1,
  ncol = 2
)
```

![](Plotting_functions_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

You can plot the mutation burden of this patient and compare it to the
rest of the cohort. This allows you to see how much this patient aligns
to the rest of the cohort in terms of clonal and subclonal mutations,
number of drivers and number of clones with drivers (summary statistics
that you obtain with).

``` r
# The data used for this plot is obtained from this function
Stats_cohort(cohort, 'CRUK0001')
```

    ## # A tibble: 1 x 7
    ##   patientID numBiopsies numMutations numDriverMutati… numClonesWithDr…
    ##   <chr>           <int>        <int>            <int>            <int>
    ## 1 CRUK0001            3         2100                7                4
    ## # … with 2 more variables: numTruncalMutations <int>,
    ## #   numSubclonalMutations <int>

``` r
# The actual scatterplot
plot_data_mutation_burden(cohort, 'CRUK0001')
```

![](Plotting_functions_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

A popular visualization for the data of a patient is the so-called
“oncoprint” where one shows for every region the presence/ absence of
a mutation. In this case we show the CCF of the mutation, and its
clone-assignment in the patient. You can also create the oncoprint
visualization of all the data available for a patient. This
visualization is done via the package `pheatmap`, and is usually quite
wide.

``` r
plot_data_oncoprint(cohort, 'CRUK0001')
```

![](Plotting_functions_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

You can also generate the marginal histogram distribution of the samples
of this patient. In this case the histograms are quite ugly because this
is exome data.

``` r
plot_data_histogram(cohort, 'CRUK0001')
```

![](Plotting_functions_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

And you have a wrap-up function `plot_data` that generates all of the
above visaulizations and assemble them in a nice figure with `ggpubr`.
You can use `plot_data` to generate a fast visualization of all relevant
information concerning the data of a patient. This code is not run, but
it would dump to PDF a single-page report of the data for each one of
the cohort patients

``` r
# 
# Not run
# 

# Output PDF 
pdf("Cohort.pdf", width = 10, height = 10)

# Loop through all patients and generate the plot 
plots = lapply(cohort$patients, plot_data, x = cohort)

# Print to PDF
lapply(plots, print)

# Close device
dev.off()
```

# Plotting trees

REVOLVER trees can be plot by S3 functions of class `rev_phylo`, when
applied to an object of class `rev_phylo` (something that is built from
the `revolver_phylogeny` constructor function. The plotting function
`plot.rev_phylo` implements 3 different types of plots:

  - (default) plot a full tree, with all its nodes nodes coloured by
    driver status, and driver ids annotated in the plot;

  - (`icon = TRUE`) plots the tree in compact form, omitting the clones
    that do not harbour any driver;

  - (`information_transfer = TRUE`) plots the Information Transfer for a
    tree, which is just the graph of edges transferred by the Transfer
    Learning fit.

The deafult plotting behaviour is called when no extra parameter is used

``` r
# We plot the top-rank tree for the first 2 patients
plot_1 = plot(Phylo(cohort, 'CRUK0001', rank = 1))
plot_2 = plot(Phylo(cohort, 'CRUK0002', rank = 1))

ggarrange(
  plot_1,
  plot_2,
  nrow = 1, 
  ncol = 2)
```

![](Plotting_functions_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

The icon version of this tree is very similar, but omits nodes that do
not harbour a driver event. This visualization is convenient when one
needs to plot several plots.

``` r
# We take the 5 top-ranking trees for CRUK0001
strip = lapply(1:5, Phylo, x = cohort, p =  'CRUK0001')
  
# Then we plot and arrange as a 1x5 matrix
ggarrange(
  plotlist = lapply(strip, plot, icon = TRUE), 
  nrow = 1, 
  ncol = 5)
```

![](Plotting_functions_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

The Information Transfer plot, where driver events are coloured with the
same principle used to plot the tree, is obtained using
`information_transfer = TRUE` as a parameter of the S3 `plot` function.
Notice that since the transfer before and after the fit likely differ:
we can check this by comparing theirs plots with these fuctions.

``` r
# Inspect what is the best solution, and obtain its rank
rank = cohort$fit$fit_table %>%
  filter(patientID == 'CRUK0001') %>%
  pull(Solution)  

# Extract the trees (before and after the fits)
before_fit = Phylo(cohort, 'CRUK0001', rank = rank, data = 'trees')
after_fit = Phylo(cohort, 'CRUK0001', rank = 1, data = 'fits')

# Plot and assemble a figure
ggarrange(
  plot(before_fit, information_transfer = TRUE),
  plot(after_fit, information_transfer = TRUE),
  nrow = 1, 
  ncol = 2)
```

![](Plotting_functions_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

Complex layouts can be trivially assembled via `ggpubr`, as before. For
instance this code visualizes side-by-side a full tree, and its
information transfer for all the cohort.

``` r
# 
# This code is not run.
# 

# Function to carry out the task
myFun = function(patient)
{
    # Get the tree
  tree = Phylo(cohort, p = patient, rank = 1)

  ggpubr::ggarrange(
    plot(tree),
    plot(tree, information_transfer = TRUE),
    ncol = 2,
    nrow = 1
  )
}

# Output PDF 
pdf("Cohort_trees.pdf", width = 10, height = 6)

# Loop through all patients and generate the plot 
plots = lapply(cohort$patients, myFun)

# Print to PDF
lapply(plots, print)

# Close device
dev.off()
```

If one seeks to visualize the score available for the trees of a
patient, there is a function that plots a barplot which shows at the
same time the tree scores and the combinations of Information Transfer
for a patient.

``` r
plot_trees_scores(cohort, 'CRUK0001')
```

    ## [[1]]

![](Plotting_functions_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

As for data, there is a wrap-up function `plot_trees` that generates all
of the above plots at once. This function is like the `myFun` function
defined above, but it also includes the strip plot and the barplot with
the scores.

``` r
plot_trees(cohort, 'CRUK0001')
```
