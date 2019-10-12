
<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/caravagn/revolver.svg?branch=master)](https://travis-ci.org/caravagn/revolver)
[![Lifecycle:
maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
<!-- badges: end -->

# revolver <a href='caravagn.github.io/REVOLVER'><img src='man/figures/logo.png' align="right" height="139" /></a>

Thee `REVOLVER` package implements the statistical model described in
[Caravagna et al;
PMID: 30171232](https://www.ncbi.nlm.nih.gov/pubmed/30171232) to
determine trajectories of repeated evolution from multi-region
sequencing data of human cancers.

The package implements functions to process data from large cohorts, and
determine different types of phylogenetic trees from the input data of
each sequenced patients. The package provides also several functions to
plot the data, and determine clusters of tumours that share repeated
trajectories; this provides a method to stratify cancer patients for the
way their tumour evolve, reconciling tumour heterogeneity both between
and within patients.

`REVOLVER` is part of the `evoverse`, a package that gathers multiple R
packages to implement Cancer Evolution analyses; see more [about
evoverse](https://caravagn.github.io/evoverse).

#### Citation

If you use `REVOLVER`, please cite:

  - G. Caravagna, Y. Giarratano, D. Ramazzoti, I. Tomlinson, T.A.
    Graham, G. Sanguinetti, A. Sottoriva. *Detecting repeated cancer
    evolution from multi-region tumor sequencing data.* Nature Methods
    15, 707–714 (2018).

#### Help and support

`REVOLVER` has its own webpage at [GitHub
pages](https://caravagn.github.io/revolver/).

-----

### Installation

You can install the released version of `REVOLVER` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("caravagn/revolver")
```

**Important news.** Since October 2019 this master branch represents a
de-novo implementation of the `REVOLVER` method. The new software
follows the [tidy paradigm](https://www.tidyverse.org/), and provides a
much more clear interface to manipulate your data, plot it and analyse
it. A number of [new
functions](https://caravagn.github.io/revolver/reference/index.html) are
now documented and many old functions have been removed. The previous
package version is mirrored in the [pre\_tibbles
branch](https://github.com/caravagn/revolver/tree/pre_tibbles) branch,
which you can still install.

``` r
# install.packages("devtools")
devtools::install_github("caravagn/revolver",  repo = 'pre_tibbles')
```

-----

#### Copyright and contacts

Giulio Caravagna, PhD. *Institute of Cancer Research, London, UK*.

  - Personal webpage:
    [https://bit.ly/2kc9E6Y](https://sites.google.com/site/giuliocaravagna/),
  - Email address: <giulio.caravagna@icr.ac.uk> and
    <gcaravagn@gmail.com>
  - Twitter feed: [@gcaravagna](https://twitter.com/gcaravagna)
  - GitHub space: [caravagn](https://github.com/caravagn)
