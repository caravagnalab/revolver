
####### Functions to check if the objetc has what it should
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

