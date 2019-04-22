
#' Title
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
CCF_parser = function(x)
{
  tk = strsplit(x, ';')[[1]]
  tk = unlist(strsplit(tk, ':'))
  
  samples = tk[seq(1, length(tk), 2)]
  
  values = tk[seq(2, length(tk), 2)]
  names(values) = samples
  
  return(values)
}
