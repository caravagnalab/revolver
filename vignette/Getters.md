Getter functions
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

In this vignette we make use of one of the cohort objects released with
the tool.

``` r
# Data from the TRACERx cohort 
load("Cohort.RData") 

# We can use S3 object functions to retrieve simple information about the plot.
# The `print` functions runs also the `revolver_check_cohort` function which 
# tells us that some patient have only 1 clone with drivers, and therefore they
# can just be expanded.
cohort
```

    ##  [ REVOLVER - Repeated Evolution in Cancer ] 
    ## 
    ##  Dataset : My REVOLVER dataset 
    ##  Cohort  : 99 patients, 65421 variants and 79 driver events. 
    ## 
    ##  Trees per patient    : YES 
    ##  Fit via TL           : YES 
    ##  REVOLVER clustering  : NO 
    ##  Jackknife statistics : NO 
    ##   
    ## For summary statistics see `?Stats_*(x)` with * = {cohort, drivers, trees, fits, clusters, ...} 
    ## 
    ## 
    ##  WARNING  Some patients have only one clone with drivers, and therefore they will just be expanded. 
    ## # A tibble: 54 x 7
    ##    patientID numBiopsies numMutations numDriverMutati… numClonesWithDr…
    ##    <chr>           <int>        <int>            <int>            <int>
    ##  1 CRUK0007            2          109                3                1
    ##  2 CRUK0010            2          141                3                1
    ##  3 CRUK0012            2          163                1                1
    ##  4 CRUK0018            4          653                4                1
    ##  5 CRUK0019            2          110                1                1
    ##  6 CRUK0021            2          221                4                1
    ##  7 CRUK0025            3          532                3                1
    ##  8 CRUK0026            2          157                4                1
    ##  9 CRUK0028            2           74                2                1
    ## 10 CRUK0029            6          399                4                1
    ## # … with 44 more rows, and 2 more variables: numTruncalMutations <int>,
    ## #   numSubclonalMutations <int>

# Access patient-level data

**Note** Since the new release of REVOLVER, we have implemented the
internal structure of the objects using the tidy data.frame
representations from [tidyverse](https://www.tidyverse.org/). Most
functions now return `tibble` data.frames that can be processed with the
`dplyr` jargon.

### Getters

We have made available several types of getters to perform queries on
the data. Getter functions for the data have a common parametrization;
for instance getter function`Drivers` takes as input

  - `x` a REVOLVER cohort object;
  - `patients` a list of patients IDs that will be used to subset the
    outputs (all by default);

<!-- end list -->

``` r
# Access all data for a patient
Data(cohort, 'CRUK0001')
```

    ## # A tibble: 2,100 x 12
    ##    id    Misc  patientID variantID is.driver is.clonal cluster cluster_size
    ##    <chr> <chr> <chr>     <chr>     <lgl>     <lgl>     <chr>          <int>
    ##  1 __mu… CRUK… CRUK0001  PKD1L1    FALSE     FALSE     1               1136
    ##  2 __mu… CRUK… CRUK0001  PCNT      FALSE     FALSE     1               1136
    ##  3 __mu… CRUK… CRUK0001  CSF2RB    FALSE     FALSE     1               1136
    ##  4 __mu… CRUK… CRUK0001  KAZALD1   FALSE     FALSE     1               1136
    ##  5 __mu… CRUK… CRUK0001  CACNA2D1  FALSE     FALSE     1               1136
    ##  6 __mu… CRUK… CRUK0001  IPO9      FALSE     FALSE     1               1136
    ##  7 __mu… CRUK… CRUK0001  PDCD10    FALSE     FALSE     1               1136
    ##  8 __mu… CRUK… CRUK0001  CALD1     FALSE     FALSE     1               1136
    ##  9 __mu… CRUK… CRUK0001  CANT1     FALSE     FALSE     1               1136
    ## 10 __mu… CRUK… CRUK0001  IFT140    FALSE     FALSE     1               1136
    ## # … with 2,090 more rows, and 4 more variables: CCF <chr>, R1 <dbl>,
    ## #   R2 <dbl>, R3 <dbl>

``` r
# Access only the drivers for a patient
Drivers(cohort, 'CRUK0001')
```

    ## # A tibble: 7 x 12
    ##   id    Misc  patientID variantID is.driver is.clonal cluster cluster_size
    ##   <chr> <chr> <chr>     <chr>     <lgl>     <lgl>     <chr>          <int>
    ## 1 __mu… CRUK… CRUK0001  NF1       TRUE      FALSE     1               1136
    ## 2 __mu… CRUK… CRUK0001  ARHGAP35  TRUE      FALSE     2                199
    ## 3 __mu… CRUK… CRUK0001  TP53      TRUE      TRUE      3                208
    ## 4 __mu… CRUK… CRUK0001  MGA       TRUE      TRUE      3                208
    ## 5 __mu… CRUK… CRUK0001  WRN       TRUE      TRUE      3                208
    ## 6 __mu… Anno… CRUK0001  EGFR      TRUE      TRUE      3                208
    ## 7 __mu… CRUK… CRUK0001  PASK      TRUE      FALSE     5                129
    ## # … with 4 more variables: CCF <chr>, R1 <dbl>, R2 <dbl>, R3 <dbl>

``` r
# Access the names of the samples for a patient
Samples(cohort, 'CRUK0001')
```

    ## [1] "R1" "R2" "R3"

``` r
# Get the list of truncal (i.e., clonal) mutations in a patient
Truncal(cohort, 'CRUK0001')
```

    ## # A tibble: 208 x 12
    ##    id    Misc  patientID variantID is.driver is.clonal cluster cluster_size
    ##    <chr> <chr> <chr>     <chr>     <lgl>     <lgl>     <chr>          <int>
    ##  1 __mu… CRUK… CRUK0001  CDH26     FALSE     TRUE      3                208
    ##  2 __mu… CRUK… CRUK0001  MYT1      FALSE     TRUE      3                208
    ##  3 __mu… CRUK… CRUK0001  STYXL1    FALSE     TRUE      3                208
    ##  4 __mu… CRUK… CRUK0001  ATXN7     FALSE     TRUE      3                208
    ##  5 __mu… CRUK… CRUK0001  LAMA5     FALSE     TRUE      3                208
    ##  6 __mu… CRUK… CRUK0001  NINL      FALSE     TRUE      3                208
    ##  7 __mu… CRUK… CRUK0001  TRY2P     FALSE     TRUE      3                208
    ##  8 __mu… CRUK… CRUK0001  GPR141    FALSE     TRUE      3                208
    ##  9 __mu… CRUK… CRUK0001  CLDN23    FALSE     TRUE      3                208
    ## 10 __mu… CRUK… CRUK0001  AMPH      FALSE     TRUE      3                208
    ## # … with 198 more rows, and 4 more variables: CCF <chr>, R1 <dbl>,
    ## #   R2 <dbl>, R3 <dbl>

``` r
# Get the list of subclonal mutations in a patient
Subclonal(cohort, 'CRUK0001')
```

    ## # A tibble: 1,892 x 12
    ##    id    Misc  patientID variantID is.driver is.clonal cluster cluster_size
    ##    <chr> <chr> <chr>     <chr>     <lgl>     <lgl>     <chr>          <int>
    ##  1 __mu… CRUK… CRUK0001  PKD1L1    FALSE     FALSE     1               1136
    ##  2 __mu… CRUK… CRUK0001  PCNT      FALSE     FALSE     1               1136
    ##  3 __mu… CRUK… CRUK0001  CSF2RB    FALSE     FALSE     1               1136
    ##  4 __mu… CRUK… CRUK0001  KAZALD1   FALSE     FALSE     1               1136
    ##  5 __mu… CRUK… CRUK0001  CACNA2D1  FALSE     FALSE     1               1136
    ##  6 __mu… CRUK… CRUK0001  IPO9      FALSE     FALSE     1               1136
    ##  7 __mu… CRUK… CRUK0001  PDCD10    FALSE     FALSE     1               1136
    ##  8 __mu… CRUK… CRUK0001  CALD1     FALSE     FALSE     1               1136
    ##  9 __mu… CRUK… CRUK0001  CANT1     FALSE     FALSE     1               1136
    ## 10 __mu… CRUK… CRUK0001  IFT140    FALSE     FALSE     1               1136
    ## # … with 1,882 more rows, and 4 more variables: CCF <chr>, R1 <dbl>,
    ## #   R2 <dbl>, R3 <dbl>

``` r
# Return the CCF entry for all the mutations of a patient,
CCF(cohort, 'CRUK0001')
```

    ## # A tibble: 2,100 x 8
    ##    id          variantID is.driver is.clonal cluster    R1    R2    R3
    ##    <chr>       <chr>     <lgl>     <lgl>     <chr>   <dbl> <dbl> <dbl>
    ##  1 __mut_id_1  PKD1L1    FALSE     FALSE     1        0.86     0     0
    ##  2 __mut_id_2  PCNT      FALSE     FALSE     1        0.86     0     0
    ##  3 __mut_id_3  CSF2RB    FALSE     FALSE     1        0.86     0     0
    ##  4 __mut_id_4  KAZALD1   FALSE     FALSE     1        0.86     0     0
    ##  5 __mut_id_5  CACNA2D1  FALSE     FALSE     1        0.86     0     0
    ##  6 __mut_id_6  IPO9      FALSE     FALSE     1        0.86     0     0
    ##  7 __mut_id_7  PDCD10    FALSE     FALSE     1        0.86     0     0
    ##  8 __mut_id_8  CALD1     FALSE     FALSE     1        0.86     0     0
    ##  9 __mut_id_9  CANT1     FALSE     FALSE     1        0.86     0     0
    ## 10 __mut_id_10 IFT140    FALSE     FALSE     1        0.86     0     0
    ## # … with 2,090 more rows

``` r
# Return the CCF entry for all the clones of a patient, the overall CCF
# values are obtained by REVOLVER from the average of CCF values across clones.
CCF_clusters(cohort, 'CRUK0001')
```

    ## # A tibble: 11 x 7
    ##    cluster nMuts is.driver is.clonal    R1    R2    R3
    ##    <chr>   <int> <lgl>     <lgl>     <dbl> <dbl> <dbl>
    ##  1 1        1136 TRUE      FALSE     0.86  0      0   
    ##  2 3         208 TRUE      TRUE      0.99  0.99   1   
    ##  3 2         199 TRUE      FALSE     0.22  0      0.95
    ##  4 4         178 FALSE     FALSE     0.89  0.13   0   
    ##  5 5         129 TRUE      FALSE     0.74  0      0.97
    ##  6 6         103 FALSE     FALSE     0     0      0.95
    ##  7 7          84 FALSE     FALSE     0     0.84   0   
    ##  8 8          22 FALSE     FALSE     0.9   0      0.14
    ##  9 9          15 FALSE     FALSE     0.3   0.14   0.95
    ## 10 10         14 FALSE     FALSE     0.83  0.755  0   
    ## 11 11         12 FALSE     FALSE     0.325 0.83   0.96

### Summary statistics for patients’ data

You can get a broad set of summary statistics about the `cohort` for a
custom set of patients. The statistcs that are available in summarised
format are patient-level (mutatational burdern, drivers etc.), and
driver-level (frequency, clonality
etc.).

``` r
# This returns patient-level statistics like the number of biopsies, overall mutations, drivers,
# clones with drivers, truncal and subclonal mutations.
# 
# This is also synonim to `Stats(cohort)`
Stats_cohort(cohort)
```

    ## # A tibble: 99 x 7
    ##    patientID numBiopsies numMutations numDriverMutati… numClonesWithDr…
    ##    <chr>           <int>        <int>            <int>            <int>
    ##  1 CRUK0001            3         2100                7                4
    ##  2 CRUK0002            3          280                7                4
    ##  3 CRUK0003            5          330                4                2
    ##  4 CRUK0004            4          323                4                2
    ##  5 CRUK0005            4         1582                6                2
    ##  6 CRUK0006            2         1856                6                3
    ##  7 CRUK0007            2          109                3                1
    ##  8 CRUK0008            2          365                6                2
    ##  9 CRUK0009            4          487                7                2
    ## 10 CRUK0010            2          141                3                1
    ## # … with 89 more rows, and 2 more variables: numTruncalMutations <int>,
    ## #   numSubclonalMutations <int>

``` r
# This returns driver-level statistics like the number of times the driver is clonal,
# subclonal, or found in general, and for quantity normalized by cohort size (i.e., the percentage)
Stats_drivers(cohort)
```

    ## # A tibble: 79 x 7
    ##    variantID numClonal p_clonal numSubclonal p_subclonal N_tot  p_tot
    ##    <chr>         <dbl>    <dbl>        <dbl>       <dbl> <dbl>  <dbl>
    ##  1 TP53             53   0.535             3      0.0303    56 0.566 
    ##  2 KRAS             24   0.242             4      0.0404    28 0.283 
    ##  3 EGFR             21   0.212             1      0.0101    22 0.222 
    ##  4 PIK3CA           20   0.202             1      0.0101    21 0.212 
    ##  5 CDKN2A           14   0.141             0      0         14 0.141 
    ##  6 SOX2             14   0.141             0      0         14 0.141 
    ##  7 KEAP1            12   0.121             0      0         12 0.121 
    ##  8 TERT             11   0.111             2      0.0202    13 0.131 
    ##  9 FGFR1             9   0.0909            0      0          9 0.0909
    ## 10 STK11             8   0.0808            0      0          8 0.0808
    ## # … with 69 more rows

The list of all patients in the cohort is accessible as
`cohort$patients`, and these functions can be run on a smaller subset of
patients.

``` r
Stats_cohort(cohort, patients = cohort$patients[1:5])
```

    ## # A tibble: 5 x 7
    ##   patientID numBiopsies numMutations numDriverMutati… numClonesWithDr…
    ##   <chr>           <int>        <int>            <int>            <int>
    ## 1 CRUK0001            3         2100                7                4
    ## 2 CRUK0002            3          280                7                4
    ## 3 CRUK0003            5          330                4                2
    ## 4 CRUK0004            4          323                4                2
    ## 5 CRUK0005            4         1582                6                2
    ## # … with 2 more variables: numTruncalMutations <int>,
    ## #   numSubclonalMutations <int>

# Access the trees in the cohort

Trees have getters similar to the data, and getters distinguish from
trees before and after the fit.

**Note** The tree fits might be slightly different from the trees before
the fit, because their Informatin Transfer is not expanded. Therefore
keep this in mind when comparing trees.

## Getters

You can to extract the tree of a patient, before its fit. This can be
one specific tree, or all of them at once. Trees before the fit are
indexed by their rank, which is obtained from the ordering of the tree
scores, which are obtained by the evaluated tree structure before the
fit.

These getters, for instance `Phylo`, take as parameter

  - `x` the cohort object;
  - `p` the patient identifier;
  - `rank` the rank of the tree to extract;
  - `data` to decide whether one wants the trees before the fit
    (`trees`), or the actual fit tree `fits`.

By logic, if you are asking for the fit trees (`data = 'fits'`), the
`rank` parameter is not considered (because there is only 1 tree fit by
REVOLVER).

``` r
# Access the top-rank tree for a patient
Phylo(cohort, 'CRUK0001', rank = 1)
```

    ##  [ REVOLVER tree - Ranked 1/100 ] 
    ## 
    ## # A tibble: 11 x 7
    ##    cluster nMuts is.driver is.clonal    R1    R2    R3
    ##    <chr>   <int> <lgl>     <lgl>     <dbl> <dbl> <dbl>
    ##  1 1        1136 TRUE      FALSE     0.86  0      0   
    ##  2 3         208 TRUE      TRUE      0.99  0.99   1   
    ##  3 2         199 TRUE      FALSE     0.22  0      0.95
    ##  4 4         178 FALSE     FALSE     0.89  0.13   0   
    ##  5 5         129 TRUE      FALSE     0.74  0      0.97
    ##  6 6         103 FALSE     FALSE     0     0      0.95
    ##  7 7          84 FALSE     FALSE     0     0.84   0   
    ##  8 8          22 FALSE     FALSE     0.9   0      0.14
    ##  9 9          15 FALSE     FALSE     0.3   0.14   0.95
    ## 10 10         14 FALSE     FALSE     0.83  0.755  0   
    ## 11 11         12 FALSE     FALSE     0.325 0.83   0.96
    ## 
    ##  Tree shape (drivers annotated)  
    ## 
    ##   \-GL
    ##    \-3 :: TP53, MGA, WRN, EGFR
    ##     |-5 :: PASK
    ##     | \-11
    ##     |  \-9
    ##     |   \-2 :: ARHGAP35
    ##     |    \-6
    ##     |     \-8
    ##     |-4
    ##     | \-1 :: NF1
    ##     |  \-10
    ##     \-7
    ## 
    ##  Information transfer  
    ## 
    ##    TP53 ---> NF1 
    ##    MGA ---> NF1 
    ##    WRN ---> NF1 
    ##    EGFR ---> NF1 
    ##    GL ---> TP53 
    ##    GL ---> MGA 
    ##    GL ---> WRN 
    ##    GL ---> EGFR 
    ##    PASK ---> ARHGAP35 
    ##    TP53 ---> PASK 
    ##    MGA ---> PASK 
    ##    WRN ---> PASK 
    ##    EGFR ---> PASK 
    ## 
    ##  Tree score 6.5321052975374e-05

``` r
# Access all trees for a patient. We use CRUK0002 because it has only 3 trees
Phylo(cohort, 'CRUK0002', rank = NULL)
```

    ## [[1]]
    ##  [ REVOLVER tree - Ranked 1/3 ] 
    ## 
    ## # A tibble: 7 x 7
    ##   cluster nMuts is.driver is.clonal    R1    R2    R3
    ##   <chr>   <int> <lgl>     <lgl>     <dbl> <dbl> <dbl>
    ## 1 1          72 TRUE      FALSE      0     0.92  0   
    ## 2 2          69 TRUE      TRUE       0.99  0.98  0.99
    ## 3 3          48 FALSE     FALSE      0     0     0.49
    ## 4 4          29 FALSE     FALSE      0.01  0.01  0.93
    ## 5 5          24 TRUE      FALSE      0.78  0     0   
    ## 6 6          23 TRUE      FALSE      0.98  0.03  0.98
    ## 7 7          15 FALSE     FALSE      0     0.41  0   
    ## 
    ##  Tree shape (drivers annotated)  
    ## 
    ##   \-GL
    ##    \-2 :: MET, TERT
    ##     |-1 :: RB1, IKZF1, KRAS
    ##     | \-7
    ##     \-6 :: EP300
    ##      |-4
    ##      | \-3
    ##      \-5 :: NF1
    ## 
    ##  Information transfer  
    ## 
    ##    MET ---> RB1 
    ##    MET ---> IKZF1 
    ##    MET ---> KRAS 
    ##    TERT ---> RB1 
    ##    TERT ---> IKZF1 
    ##    TERT ---> KRAS 
    ##    GL ---> MET 
    ##    GL ---> TERT 
    ##    EP300 ---> NF1 
    ##    MET ---> EP300 
    ##    TERT ---> EP300 
    ## 
    ##  Tree score 0.6 
    ## 
    ## [[2]]
    ##  [ REVOLVER tree - Ranked 2/3 ] 
    ## 
    ## # A tibble: 7 x 7
    ##   cluster nMuts is.driver is.clonal    R1    R2    R3
    ##   <chr>   <int> <lgl>     <lgl>     <dbl> <dbl> <dbl>
    ## 1 1          72 TRUE      FALSE      0     0.92  0   
    ## 2 2          69 TRUE      TRUE       0.99  0.98  0.99
    ## 3 3          48 FALSE     FALSE      0     0     0.49
    ## 4 4          29 FALSE     FALSE      0.01  0.01  0.93
    ## 5 5          24 TRUE      FALSE      0.78  0     0   
    ## 6 6          23 TRUE      FALSE      0.98  0.03  0.98
    ## 7 7          15 FALSE     FALSE      0     0.41  0   
    ## 
    ##  Tree shape (drivers annotated)  
    ## 
    ##   \-GL
    ##    \-2 :: MET, TERT
    ##     \-1 :: RB1, IKZF1, KRAS
    ##      |-6 :: EP300
    ##      | |-4
    ##      | | \-3
    ##      | \-5 :: NF1
    ##      \-7
    ## 
    ##  Information transfer  
    ## 
    ##    MET ---> RB1 
    ##    MET ---> IKZF1 
    ##    MET ---> KRAS 
    ##    TERT ---> RB1 
    ##    TERT ---> IKZF1 
    ##    TERT ---> KRAS 
    ##    GL ---> MET 
    ##    GL ---> TERT 
    ##    EP300 ---> NF1 
    ##    RB1 ---> EP300 
    ##    IKZF1 ---> EP300 
    ##    KRAS ---> EP300 
    ## 
    ##  Tree score 0.0666666666666667 
    ## 
    ## [[3]]
    ##  [ REVOLVER tree - Ranked 3/3 ] 
    ## 
    ## # A tibble: 7 x 7
    ##   cluster nMuts is.driver is.clonal    R1    R2    R3
    ##   <chr>   <int> <lgl>     <lgl>     <dbl> <dbl> <dbl>
    ## 1 1          72 TRUE      FALSE      0     0.92  0   
    ## 2 2          69 TRUE      TRUE       0.99  0.98  0.99
    ## 3 3          48 FALSE     FALSE      0     0     0.49
    ## 4 4          29 FALSE     FALSE      0.01  0.01  0.93
    ## 5 5          24 TRUE      FALSE      0.78  0     0   
    ## 6 6          23 TRUE      FALSE      0.98  0.03  0.98
    ## 7 7          15 FALSE     FALSE      0     0.41  0   
    ## 
    ##  Tree shape (drivers annotated)  
    ## 
    ##   \-GL
    ##    \-2 :: MET, TERT
    ##     \-1 :: RB1, IKZF1, KRAS
    ##      \-7
    ##       \-6 :: EP300
    ##        |-4
    ##        | \-3
    ##        \-5 :: NF1
    ## 
    ##  Information transfer  
    ## 
    ##    MET ---> RB1 
    ##    MET ---> IKZF1 
    ##    MET ---> KRAS 
    ##    TERT ---> RB1 
    ##    TERT ---> IKZF1 
    ##    TERT ---> KRAS 
    ##    GL ---> MET 
    ##    GL ---> TERT 
    ##    EP300 ---> NF1 
    ##    RB1 ---> EP300 
    ##    IKZF1 ---> EP300 
    ##    KRAS ---> EP300 
    ## 
    ##  Tree score 0.0666666666666667

Notice that in the printing of a tree to screen you can immediately see
the Information Transfer (IT) for the driver genes. In general, you can
access the IT of a tree with another getter, which takes as extra
parameter `type` in order to return either the transfer across drivers,
or across clones annotated in a tree.

``` r
# Information Transfer for the drivers, top-ranking tree
ITransfer(cohort, "CRUK0001", rank = 1, type = 'drivers')
```

    ## # A tibble: 13 x 2
    ##    from  to      
    ##    <chr> <chr>   
    ##  1 TP53  NF1     
    ##  2 MGA   NF1     
    ##  3 WRN   NF1     
    ##  4 EGFR  NF1     
    ##  5 GL    TP53    
    ##  6 GL    MGA     
    ##  7 GL    WRN     
    ##  8 GL    EGFR    
    ##  9 PASK  ARHGAP35
    ## 10 TP53  PASK    
    ## 11 MGA   PASK    
    ## 12 WRN   PASK    
    ## 13 EGFR  PASK

``` r
# Information Transfer for the clones, top-ranking tree
ITransfer(cohort, "CRUK0001", rank = 1, type = 'clones')
```

    ## # A tibble: 4 x 2
    ##   from  to   
    ##   <chr> <chr>
    ## 1 3     1    
    ## 2 GL    3    
    ## 3 5     2    
    ## 4 3     5

Fit trees can be accessed using the `data` argument. Essentially this is
like before, but does not require specifying a `rank` parameter.

``` r
# Access the fit tree for a patient
Phylo(cohort, 'CRUK0001', data = 'fits')
```

    ##  [ REVOLVER tree - Ranked 3/100 - Information Transfer expanded via Transfer Learning ] 
    ## 
    ## # A tibble: 11 x 7
    ##    cluster nMuts is.driver is.clonal    R1    R2    R3
    ##    <chr>   <int> <lgl>     <lgl>     <dbl> <dbl> <dbl>
    ##  1 1        1136 TRUE      FALSE     0.86  0      0   
    ##  2 3         208 TRUE      TRUE      0.99  0.99   1   
    ##  3 2         199 TRUE      FALSE     0.22  0      0.95
    ##  4 4         178 FALSE     FALSE     0.89  0.13   0   
    ##  5 5         129 TRUE      FALSE     0.74  0      0.97
    ##  6 6         103 FALSE     FALSE     0     0      0.95
    ##  7 7          84 FALSE     FALSE     0     0.84   0   
    ##  8 8          22 FALSE     FALSE     0.9   0      0.14
    ##  9 9          15 FALSE     FALSE     0.3   0.14   0.95
    ## 10 10         14 FALSE     FALSE     0.83  0.755  0   
    ## 11 11         12 FALSE     FALSE     0.325 0.83   0.96
    ## 
    ##  Tree shape (drivers annotated)  
    ## 
    ##   \-GL
    ##    \-3 :: TP53, MGA, WRN, EGFR
    ##     |-4
    ##     | \-1 :: NF1
    ##     |  \-10
    ##     |   \-5 :: PASK
    ##     |    \-11
    ##     |     \-9
    ##     |      \-2 :: ARHGAP35
    ##     |       \-6
    ##     |        \-8
    ##     \-7
    ## 
    ##  Information transfer  
    ## 
    ##    GL ---> EGFR 
    ##    GL ---> WRN 
    ##    GL ---> MGA 
    ##    EGFR ---> TP53 
    ##    WRN ---> TP53 
    ##    TP53 ---> NF1 
    ##    MGA ---> NF1 
    ##    NF1 ---> PASK 
    ##    PASK ---> ARHGAP35 
    ## 
    ##  Tree score 6.5321052975374e-05

``` r
# Information Transfer for the drivers, top-ranking tree. Notice that this is different
# from the result of the above call, because the transfer after fitting is expanded
ITransfer(cohort, "CRUK0001", rank = 1, type = 'drivers', data = 'fits')
```

    ## # A tibble: 9 x 2
    ##   from  to      
    ##   <chr> <chr>   
    ## 1 GL    EGFR    
    ## 2 GL    WRN     
    ## 3 GL    MGA     
    ## 4 EGFR  TP53    
    ## 5 WRN   TP53    
    ## 6 TP53  NF1     
    ## 7 MGA   NF1     
    ## 8 NF1   PASK    
    ## 9 PASK  ARHGAP35

### Summary statistics for trees and fits

There are getters for summary statistics that work for trees and fits,
with the same principles fo the getters for the data discussed
above

``` r
# This returns patient-level statistics for the trees available in a patient. The tibble reports
# whether the patient has trees annotated, the total number of trees, their minimum and maximum
# scores mutations and the total number of differnet combinations of Information Transfer for 
# the available trees.
Stats_trees(cohort)
```

    ## # A tibble: 99 x 6
    ##    patientID hasTrees numTrees  maxScore  minScore combInfTransf
    ##    <chr>     <lgl>       <int>     <dbl>     <dbl>         <int>
    ##  1 CRUK0001  TRUE          100 0.0000653 0.0000218             8
    ##  2 CRUK0002  TRUE            3 0.6       0.0667                2
    ##  3 CRUK0003  TRUE            1 0.8       0.8                   1
    ##  4 CRUK0004  TRUE           18 0.02      0.005                 1
    ##  5 CRUK0005  TRUE          100 0.000126  0.0000628             1
    ##  6 CRUK0006  TRUE          100 0.000347  0.0000868             3
    ##  7 CRUK0007  TRUE            1 1         1                     1
    ##  8 CRUK0008  TRUE            1 1         1                     1
    ##  9 CRUK0009  TRUE           55 0.00617   0.000193              1
    ## 10 CRUK0010  TRUE            1 1         1                     1
    ## # … with 89 more rows

``` r
# This returns the same table of above, but with some extended information on the fits (like the fit rank, etc)
# Stats_fits(cohort)
```

# Plotting functions

``` r
plot(Phylo(cohort, 'CRUK0001', rank = 1))
```

![](Getters_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
plot(Phylo(cohort, 'CRUK0001', rank = 1), information_transfer = TRUE)
```

![](Getters_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
plot(Phylo(cohort, 'CRUK0001', rank = 1), icon = TRUE)
```

![](Getters_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->
