

jamPDF = function(in.files, out.file = 'jamPDF.pdf', layout = '3x3',
                  delete.original = TRUE, hide.output = TRUE, crop.white = TRUE,
                  page = NULL)
{
  in.files = in.files[sapply(in.files, file.exists)]
  if(length(in.files) == 0) return()

  cmd = paste('pdfjam ',
              paste(in.files, collapse = ' '),
              ' --nup ', layout, ' --landscape --outfile ', out.file, sep = ' ')

  aa = system(cmd, intern = hide.output, ignore.stderr = TRUE)

  if(crop.white)
    aa = system(
      paste('pdfcrop --margins "3 3 3 3"', out.file, out.file),
      intern = hide.output, ignore.stderr = TRUE)

  if(!is.null(page))
  {
    page = R.utils::capitalize(page)
    aa = system(
      paste(
        paste(GIT,'/bin/pdfScale.sh', sep = ''),  '-v -r', page, out.file),
      intern = hide.output, ignore.stderr = TRUE)

    f = gsub(out.file, pattern = '.pdf', replacement = '')
    file.rename(paste(f, page, 'pdf', sep ='.'), paste(f, 'pdf', sep ='.'))
  }

  if(delete.original) file.remove(setdiff(in.files, out.file))

  invisible(return())
}


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid::grid.newpage()
    grid::pushViewport(
      grid::viewport(
        layout = grid::grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = grid::viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
    }
  }
}


scols = function(v, palette = 'Set1', alpha = 1)
{
  pmax = RColorBrewer::brewer.pal.info[palette, 'maxcolors']
  colors = RColorBrewer::brewer.pal(n = pmax, palette)

  if(length(v) > pmax) colors = colorRampPalette(colors)(length(v))

  colors = colors[1:length(v)]

  add.alpha = Vectorize(add.alpha, vectorize.args = 'col')
  add.alpha(colors, alpha)

  names(colors) = v

  return(colors)
}

add.alpha =  function(col, alpha=1) {
  apply(sapply(col, col2rgb)/255, 2,
        function(x)
          rgb(x[1], x[2], x[3], alpha=alpha))
}



dodev = function(width = 20, height = 15, device = 'quartz'){

  if(device == 'quartz') quartz(width = width, height = height)
  else dev.new(width = width, height = height)
}

udodev = function(file = NA){
  if(!is.na(file)) dev.copy2pdf(file = file)
  dev.off()
}


prtTit = function(title) {
  ncharT = nchar(title)

  line = cyan(paste(rep('-', ncharT), collapse = ''))

  cat(paste('\n', line, '\n', sep = ''))
  cat(cyan(title))
  cat(paste('\n', line, '\n', sep = ''))
}

prtHdr = function(header, ..., format = '') {
  cat(bgBlue(paste(" [", header, "] \n")))

  toPrint = list(...)
  ntoPrint = length(toPrint)/2

  for(l in seq(1, length(toPrint), 2)) cat(format, cyan(toPrint[[l]]), toPrint[[l+1]], '\n')
}


