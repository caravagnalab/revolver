#' Palette of few distinct colours.
#' 
#' @description ColorRamp palette extrapolating \code{Set1} palette of \code{RColorBrewer} package. The original \code{Set1} palette is
#' sampled for 9 colours.
#' 
#' @return A function that applied to a number returns the required colours.
#' @export
#'
#' @examples
#' distinct_palette_few(2)
distinct_palette_few = function()
{
  colorRampPalette(RColorBrewer::brewer.pal(n = 9, "Set1"))
}

#' Palette of \code{n} distinct colours.
#' 
#' @description ColorRamp palette extrapolating all qualitative (\code{qual}) palettes from
#' the \code{RColorBrewer} package. This function requires the number of required colours
#' as input parameter.
#' 
#' @param n The number of required colours.
#'
#' @return A function that applied to a number returns the required colours.
#' @export
#'
#' @examples
#' distinct_palette_many(2)
distinct_palette_many = function(n) {
  # Qualitative colors
  qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual', ]
  # All colors
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  colorRampPalette(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))) (n)
}


#' Gradient palette.
#' 
#' @description ColorRamp palette extrapolating \code{YlGnBu} palette of
#' \code{RColorBrewer} package. The original \code{YlGnBu} palette is
#' sampled for 9 colours.
#'
#' @return A function that applied to a number returns the required colours.
#' @export
#'
#' @examples
#' gradient_palette(2)
gradient_palette = function(){
  colorRampPalette(RColorBrewer::brewer.pal(n = 9, "YlGnBu"))
}
