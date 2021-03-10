
<!-- badges: start -->

[![R-CMD-check](https://github.com/caravagnalab/revolver/workflows/R-CMD-check/badge.svg)](https://github.com/caravagnalab/revolver/actions)
[![pkgdown](https://github.com/caravagnalab/revolver/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/caravagnalab/revolver/actions/workflows/pkgdown.yaml)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-green.svg)](https://www.tidyverse.org/lifecycle/#stable)
<!-- badges: end -->

# revolver <a href='caravagn.github.io/revolver'><img src='man/figures/logo.png' align="right" height="139" /></a>

The `revolver` package implements the statistical model described in
[Caravagna et al; PMID:
30171232](https://www.ncbi.nlm.nih.gov/pubmed/30171232) to determine
trajectories of repeated evolution from multi-region sequencing data of
human cancers.

The package implements functions to process data from large cohorts, and
determine different types of phylogenetic trees from the input data of
each sequenced patients. The package provides also several functions to
plot the data, and determine clusters of tumours that share repeated
trajectories; this provides a method to stratify cancer patients for the
way their tumour evolve, reconciling tumour heterogeneity both between
and within patients.

#### Citation

If you use `revolver`, please cite:

-   G. Caravagna, Y. Giarratano, D. Ramazzoti, I. Tomlinson, T.A.
    Graham, G. Sanguinetti, A. Sottoriva. *Detecting repeated cancer
    evolution from multi-region tumor sequencing data.* Nature Methods
    15, 707–714 (2018).

[![](https://img.shields.io/badge/doi-10.1038/s41592--018--0108--x-red.svg)](https://doi.org/10.1038/s41592-018-0108-x)

#### Help and support

[![](https://img.shields.io/badge/GitHub%20Pages-https://caravagnalab.github.io/revolver/-yellow.svg)](https://caravagnalab.github.io/revolver)

------------------------------------------------------------------------

### Installation

You can install the released version of `revolver` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("caravagnalab.github.io/revolver")
```

------------------------------------------------------------------------

#### Copyright and contacts

Giulio Caravagna. Cancer Data Science (CDS) Laboratory.

[![](https://img.shields.io/badge/Email-gcaravagn@gmail.com-steelblue.svg)](mailto:gcaravagn@gmail.com)
[![](https://img.shields.io/badge/CDS%20Lab%20Github-caravagnalab-seagreen.svg)](https://github.com/caravagnalab)
[![](https://img.shields.io/badge/CDS%20Lab%20webpage-https://www.caravagnalab.org/-red.svg)](https://www.caravagnalab.org/)
