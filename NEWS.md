# revolver 1.0.1

* Fixed missing package dependencies in `DESCRIPTION`: added `cli`, `cowplot`, `reshape2`, and `broom` to `Imports`.
* Removed `pheatmap`, `R6`, and `clisymbols` from `Imports` (unused); moved `evoverse.datasets` to `Suggests`.
* Replaced deprecated `dplyr::progress_estimated()` (removed in dplyr 1.0.0) with `cli::cli_progress_bar()` in `revolver_fit`, `revolver_cluster`, `remove_drivers`, and `utils_jackknife`.
* Removed `require()` call inside `plot_patient_oncoprint` (non-standard in package code).
* Fixed typos in `DESCRIPTION` description field.
* Fixed incorrect install URL in `README.md`.
* Updated GitHub Actions workflows to current versions (`actions/checkout@v4`, `r-lib/actions/*@v2`).
* Fixed duplicate `@param x` roxygen entries in `compute_clone_trees` and `compute_mutation_trees`; added missing `@return` and `@param` documentation.
* Fixed `revolver_cluster` example referencing non-existent `revolver_evo_distance()`.
* Added `^old_vignette$` to `.Rbuildignore`.
* Removed large commented-out dead code block from `revolver_jackknife.R`.
* Updated `Date` to 2026-06-29.

# revolver 1.0.0

## revolver 0.3.0 

* Re-coding of all the internal structure of the cohorts using tibbles and in general the tidy data paradigm;
* New interface to access the data, based on ad hoc getter functions;
* All plotting functions have been re-implemented using ggplot to gain more flexibility, and remove the PDF jamming routines;
* Clone trees from CCF data have been extracted and included in a new dedicated package [caravagn.github.io/ctree](caravagn.github.io/ctree);
* Added a `NEWS.md` file to track changes to the package.

## revolver 0.2.0 

* Bugs fixing, and removal of the actual code taken from CloneEvol to increase clarity and improve performance.
* jamPDF is now associated to an option (default = FALSE) to increase compatibility with MS Windows systems.

## revolver 0.1.0 

* Released the first version of the tool, together with the paper PMID: 30171232.
