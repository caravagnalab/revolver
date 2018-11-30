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
#' @import grDevices
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
  toJam = getOption("revolver.jamPDF", default = FALSE)
  if(!toJam)
  {
    message("jamPDF is disabled by default. Use \n\n\toptions('revolver.jamPDF' = TRUE)\n\nto use this funcitonality and merge PDFs scripts; jamPDF requires 'pdfjam' installed to work.")
    return(invisible(1))
  }

  # DISABLED
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


