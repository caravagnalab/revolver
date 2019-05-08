# Distinct colors if less than 9
distinct_palette_few = colorRampPalette(RColorBrewer::brewer.pal(n = 9, "Set1"))

# Distinct colors with loads of colors
distinct_palette_many = function(n) {
  # Qualitative colors
  qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual', ]
  # All colors
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  colorRampPalette(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))) (n)
}


# .....
gradient_palette = colorRampPalette(RColorBrewer::brewer.pal(n = 9, "YlGnBu"))
