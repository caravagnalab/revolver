#' Palette of few distinct colours.
#' 
#' @description ColorRamp palette extrapolating \code{Set1} palette of \code{RColorBrewer} package. The original \code{Set1} palette is
#' sampled for 9 colours. This function can have the number of required colours
#' as input parameter; if this does not happen a function is returned.
#' 
#' @param n The number of required colours. By default this is NULL.
#' @family Plotting functions
#' @return Either a function that applied to a number returns the required colours,
#' or the actual sampling function.
#' @export
#'
#' @examples
#' distinct_palette_few(2)
distinct_palette_few = function(n = NULL)
{
  f = colorRampPalette(RColorBrewer::brewer.pal(n = 9, "Set1"))
  
  if(is.null(n)) return(f)
  else return(f(n))
}

#' Palette of \code{n} distinct colours.
#' 
#' @description ColorRamp palette extrapolating all qualitative (\code{qual}) palettes from
#' the \code{RColorBrewer} package. This function can have the number of required colours
#' as input parameter; if this does not happen a function is returned.
#' 
#' @param n The number of required colours. By default this is NULL.
#'
#' @family Plotting functions
#' @return Either a function that applied to a number returns the required colours,
#' or the actual sampling function.
#' 
#' @export
#'
#' @examples
#' distinct_palette_many(2)
distinct_palette_many = function(n = NULL) {
  # Qualitative colors
  qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual', ]
  # All colors
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  f = colorRampPalette(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
  
  if(is.null(n)) return(f)
  else return(f(n))
}


#' Gradient palette.
#' 
#' @description ColorRamp palette extrapolating \code{YlGnBu} palette of
#' \code{RColorBrewer} package. The original \code{YlGnBu} palette is
#' sampled for 9 colours.
#' 
#' @param n The number of required colours. By default this is NULL.
#'
#' @family Plotting functions
#' @return Either a function that applied to a number returns the required colours,
#' or the actual sampling function.
#' 
#' @export
#'
#' @examples
#' gradient_palette(2)
gradient_palette = function(n = NULL) {
  
  f = colorRampPalette(RColorBrewer::brewer.pal(n = 9, "YlGnBu"))
  
  if(is.null(n)) return(f)
  else return(f(n))
}
