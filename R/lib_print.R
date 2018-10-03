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
multiplot <-
  function(...,
           plotlist = NULL,
           file,
           cols = 1,
           layout = NULL) {

    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)

    numPlots = length(plots)

    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
      # Make the panel
      # ncol: Number of columns of plots
      # nrow: Number of rows needed, calculated from # of cols
      layout <- matrix(seq(1, cols * ceiling(numPlots / cols)),
                       ncol = cols,
                       nrow = ceiling(numPlots / cols))
    }

    if (numPlots == 1) {
      print(plots[[1]])

    } else {
      # Set up the page
      grid::grid.newpage()
      grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow(layout), ncol(layout))))

      # Make each plot, in the correct location
      for (i in 1:numPlots) {
        # Get the i,j matrix positions of the regions that contain this subplot
        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

        print(
          plots[[i]],
          vp = grid::viewport(
            layout.pos.row = matchidx$row,
            layout.pos.col = matchidx$col
          )
        )
      }
    }
  }


scols = function(v, palette = 'Set1', alpha = 1)
{
  pmax = RColorBrewer::brewer.pal.info[palette, 'maxcolors']
  colors = RColorBrewer::brewer.pal(n = pmax, palette)

  if (length(v) > pmax)
    colors = colorRampPalette(colors)(length(v))

  colors = colors[1:length(v)]

  add.alpha = Vectorize(add.alpha, vectorize.args = 'col')
  add.alpha(colors, alpha)

  names(colors) = v

  return(colors)
}

add.alpha =  function(col, alpha = 1) {
  apply(sapply(col, col2rgb) / 255, 2,
        function(x)
          rgb(x[1], x[2], x[3], alpha = alpha))
}



dodev = function(width = 20,
                 height = 15,
                 device = 'quartz') {
  if (device == 'quartz')
    quartz(width = width, height = height)
  else
    dev.new(width = width, height = height)
}

udodev = function(file = NA) {
  if (!is.na(file))
    dev.copy2pdf(file = file)
  dev.off()
}


mylayout.on = function(file = NA,
                       nplots = 1,
                       size = c(10, 10),
                       cex = 1)
{
  cur.params = par()
  cur.size = dev.size()

  if (!is.na(file))
    pdf(file, width = size[1] * cex, height = size[2] * cex)
  else
  {
    sq = ceiling(sqrt(nplots))

    # size = size * sq
    # if(any(size > cur.size))
    # {
    #   dev.off()
    #
    #   dev.new(noRStudioGD = TRUE,
    #           width = size[1],
    #           height = size[2])
    # }

    if (nplots > 1)
      par(mfrow = c(sq, sq))
  }

  dev = list(
    cur.params = cur.params,
    cur.size = cur.size,
    file = file,
    nplots = nplots
  )

  return(dev)
}


mylayout.off = function(config) {
  if (!is.na(config$file))
    dev.off()
  else {
    if (config$nplots > 1)
      par(mfrow = config$cur.params$mfrow)
    #
    #     dev.off()
    #
    #     dev.new(noRStudioGD = TRUE,
    #             width = config$cur.size[1],
    #             height = config$cur.size[2])
  }
}
