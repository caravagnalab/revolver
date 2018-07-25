


#' PDF jamming function
#'
#' @description
#' This function allows to combine PDFs in a simple way. It uses pdfjam to assemble files,
#' and pdfScale.sh to resize the output. The script is included in inst/bin, but pdfjam no.
#'
#' @param in.files A bunch of input PDF file names
#' @param out.file Output file, cannot be NA
#' @param layout A layout to assemble PDFs. For isntance "1x1" is one PDF per page, etc.
#' @param delete.original TRUE to delete the input files
#' @param hide.output TRUE to avoid showing to screen the output of system() calls
#' @param crop.white TRUE to trim each margin by a factor 3 (white space removal)
#' @param page Set it to any like "a4" or "a3" etc to resize each PDF to that format. Leave it to "none"
#' to avoid this step
#' @param ignore.stderr TRUE to ignore standard error
#'
#' @importFrom grDevices col2rgb colorRampPalette dev.copy2pdf dev.new
#' @importFrom grDevices dev.off dev.size pdf quartz rgb
#' 
#' @return None
#' @export
#'
jamPDF = function(in.files,
                  out.file = 'jamPDF.pdf',
                  layout = '3x3',
                  delete.original = TRUE,
                  crop.white = TRUE,
                  page = 'none',
                  hide.output = TRUE,
                  ignore.stderr = TRUE)
{
  in.files = in.files[sapply(in.files, file.exists)]

  if(length(in.files) == 0) stop("All the input files that you asked to jam are missing -- check input!")

  pio::pioHdr(paste("REVOLVER jamPDF to", out.file),
              toPrint = c(
                `Input files` = paste(in.files, collapse = ', '),
                `PDF layout` = layout,
                `PDF page type` = page,
                `PDF crop white margins` = crop.white,
                `PDF delete input files` = delete.original),
              prefix = '\t -'
              )

  stopifnot(!is.na(out.file))


  cmd = paste(
    'pdfjam ',
    paste(in.files, collapse = ' '),
    ' --nup ',
    layout,
    ' --landscape --outfile ',
    out.file,
    sep = ' '
  )


  # pio::pioTit("Assembling PDFs")
  aa = system(cmd, intern = hide.output, ignore.stderr = ignore.stderr)

  if (crop.white)
  {
    # pio::pioTit("Cropping white margins (3 3 3 3)")

    aa = system(
      paste('pdfcrop --margins "3 3 3 3"', out.file, out.file),
      intern = hide.output,
      ignore.stderr = ignore.stderr
    )
  }


  if (page != 'none')
  {
    page = R.utils::capitalize(page)

    # pio::pioTit(paste("Resizing pages to", page))


    installation = find.package('revolver')
    cmd = paste(paste0(installation, '/bin/pdfScale.sh'),
                '-v -r',
                page,
                out.file)


    aa = system(cmd, intern = hide.output, ignore.stderr = ignore.stderr)


    f = gsub(out.file, pattern = '.pdf', replacement = '')
    file.rename(paste(f, page, 'pdf', sep = '.'), paste(f, 'pdf', sep =
                                                          '.'))
  }

  if (delete.original)
    file.remove(setdiff(in.files, out.file))

  invisible(NULL)
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
