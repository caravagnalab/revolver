
####### Functions to check if the objetc has what it should

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# To test if the object has some consistency internally
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

stop_not_revolver_object = function(x)
{
  pass = any(sapply(c('rev_cohort', 'rev_cohort_fit'), inherits, x = x))

  if(!pass)
    stop("Input object is not a REVOLVER cohort, aborting.")
}

stop_invalid_patient = function(x, p)
{
  pass = p %in% x$patients

  if(!pass)
    stop("Patient", p, "does not exist in this REVOLVER cohort, aborting.")
}

has_patient_trees = function(x, p = NULL, rank = NULL)
{
  if (is.null(p))
    return(!is.null(x$phylogenies))

  if (is.null(rank))
    return(p %in% names(x$phylogenies))
  else
    return(rank <= length(x$phylogenies[[p]]))
}

has_fits = function(x, p = NULL)
{
  if (is.null(p))
    return(!is.null(x$fit))

  models = x$fit$phylogenies[p]

  return(!all(sapply(models, is.null)))
}
# has_fits = Vectorize(has_fits, vectorize.args = 'p', SIMPLIFY = TRUE)

has_clusters= function(x, p = NULL)
{
  if (is.null(p))
    return(!is.null(x$cluster))

  return(p %in% names(x$cluster$fits$labels))
}

# Simpler tests
obj_has_clusters = function(x)
{
  if(is.null(x$cluster))
    stop('Cannot proceed without having computed clusters -- are you calling the right function?')
}

obj_has_fit = function(x)
{
  if(is.null(x$fit))
    stop('Cannot proceed without having computed fit models -- are you calling the right function?')
}

obj_has_evodistance = function(x)
{
  if(is.null(x$cluster$distances))
    stop('Cannot proceed without having computed the evolutionary distance -- are you calling the right function?')
}

obj_has_jackknife = function(x)
{
  if(is.null(x$jackknife))
    stop('Cannot proceed without having computed the jackknife -- are you calling the right function?')
}

