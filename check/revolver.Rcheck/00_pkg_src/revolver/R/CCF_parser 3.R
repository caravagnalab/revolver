
#' Parser function for the  required \code{CCF} input field.
#' 
#' @description This function can parse a string
#' in the format \code{R1:0.1;R2:0.4;R3:0.9} to extract
#' a 3-dimensional named vector with names \code{R1}, \code{R2} and
#' \code{R3}, and values \code{0.1}, \code{0.4} and \code{0.9}. 
#' 
#' This function can be used when parsing an input dataset in 
#' \code{revolver_cohort}, the function that creates a \code{REVOLVER} cohort. 
#' 
#' @note This function will not convert the types of the parsed entries.
#' 
#' @param x The string to parse, the accepted format for this function
#' is something like \code{R1:0.1;R2:0.4;R3:0.9}.
#'
#' @family Cohort creation
#'
#' @return A parsed named vector.
#' 
#' @export
#'
#' @examples
#' # The output is a string 3-dimensional vector  
#' CCF_parser("R1:0.1;R2:0.4;R3:0.9")
CCF_parser = function(x)
{
  tk = strsplit(x, ';')[[1]]
  tk = unlist(strsplit(tk, ':'))
  
  samples = tk[seq(1, length(tk), 2)]
  
  values = tk[seq(2, length(tk), 2)]
  names(values) = samples
  
  return(values)
}
