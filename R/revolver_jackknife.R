#' Compute clusters stability via the jackknife.
#'
#' @description
#'
#' For a set of clusters computed via \code{\link{revolver_cluster}}, you can compute
#' their stability via a jackknife. routine. This funcion runs a kind of bootstrap
#' routine where subset of patients - a desired number - is removed from the cohort and
#' before re-computing the clusters. In this way, the co-clustering probability of each
#' patient is computed, which leads to a mean clustering stability for each one of the
#' original set of clusters and a frequency for the inference of a particular evolutionary
#' trajectory.
#'
#' A number of functions are available to plot the results from this jackknife analysis.
#' Note that in general if you require a large number of runs (i.e., resamples), this
#' computation can take some time. This implementation leverages on the \code{easypar}
#' package to run in parallel all the re-runs, therefore we suggest to run it on a
#' multi-core machine to appreciate a speed up in the computations.
#'
#' @param x A cohort object where fit and clusters have been computed.
#' @param resamples Number of jackknife samples.
#' @param leave.out A number in \code{(0,1)} for the proportion of patients to leave out in each jackknife iteration.
#' @param options.fit List of parameters for fitting models. See  \code{\link{revolver_fit}}.
#' @param options.clustering List of parameters for clustering. See \code{\link{revolver_cluster}}.
#' @param ... Additional parameters forwarded to \code{\link[easypar]{run}}.
#'
#' @return A cohort where a new jackknife field contains result from this analysis
#' @export
#' @import easypar
#'
#' @importFrom utils combn
#'
#' @examples
#' \dontrun{
#' data('TRACERx_NEJM_2017_REVOLVER', package = 'evoverse.datasets')
#' x = revolver_jackknife(TRACERx_NEJM_2017_REVOLVER)
#' }
#'
revolver_jackknife = function(x,
                              resamples = 100,
                              leave.out = 0.1,
                              options.fit = list(initial.solution = NA, max.iterations = 10, n = 10),
                              options.clustering = list(min.group.size = 3, hc.method = 'ward', split.method = 'cutreeHybrid'),
                              ...
)
{
  obj_has_fit(x)
  obj_has_clusters(x)

  stopifnot(leave.out > 0 & leave.out < 1)
  stopifnot(resamples >= 1)

  # warning("October 2019 - this is to be implemented with the new easypar.")
  # return(cohort)

  patients = x$patients
  npatients = length(patients)

  pio::pioHdr('REVOLVER Jackknife')

  # Inputs samples lists
  jf_cohorts = lapply(1:resamples, function(w) {
    n = npatients * (1 - leave.out)
    sample(patients, n)
    })
  
  # Fitting function
  fitting_function = function(i) {
    y = x
    
    # ellipsis parameters
    options_fit = unpack(..., params = c('initial.solution', 'max.iterations', 'n'), fun = revolver_fit)
    options_clustering = unpack(..., params = c('hc.method', 'split.method', 'min.group.size'), fun = revolver_cluster)
    
    # patients to remove - this also remove non recurring driver events
    y = remove_patients(y, setdiff(x$patients, jf_cohorts[[i]]))

    # fitting with the required parameters, without parallel and progress bar
    y_fit = revolver_fit(y,
                         initial.solution = options_fit$initial.solution,
                         max.iterations = options_fit$max.iterations,
                         n = options_fit$n,
                         parallel = FALSE, 
                         progress_bar = FALSE)

    # fitting with the required parameters
    y_fit_clustering = revolver_cluster(y_fit,
                                        hc.method = options_clustering$hc.method,
                                        split.method = options_clustering$split.method,
                                        min.group.size = options_clustering$min.group.size)

    return(y_fit_clustering)
  }

  # jackknife resamples
  jackknife_resamples = easypar::run(
    FUN = fitting_function,
    PARAMS = lapply(1:resamples, list),
    export = ls(globalenv(), all.names = TRUE),
    packages = 'revolver',
    parallel = TRUE
  )

  # Polish errors if any
  nerrs = easypar::numErrors(jackknife_resamples)
  if(nerrs > 0) warning("Errors in the jackknife procedure, returning less results than required...")

  jackknife_resamples = easypar::filterErrors(jackknife_resamples)

  # Add the resamples and other data
  x$jackknife$results = jackknife_resamples
  x$jackknife$params = list(resamples = resamples, leave.out = leave.out)

  x = analyse_jackknife(x)

  return(x)
}

unpack = function(..., params, fun) 
{
  parameters = list(...)
  
  # values with default
  values = formals(fun)
  values = values[sapply(values, function(x) x != '')]
  
  # Specified parameters
  w_params = which(params %in% parameters)
  parameters = parameters[params] %>% as.vector
  
  # Substitute the default with the one specified
  w_specified = which(names(values) %in% names(parameters))
  values[w_specified] = parameters[names(values)[w_specified]]
  
  values
}

