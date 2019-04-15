revolver_version = function()
{
  REVOLVER_VERSION_NUMBER = "1.0"
  REVOLVER_VERSION_CODNAME = "\"Haggis and tatties\""

  REVOLVER_VERSION = paste(REVOLVER_VERSION_NUMBER, REVOLVER_VERSION_CODNAME, sep = ' ~ ')

  REVOLVER_AUTHOR = "\"Giulio Caravagna (ICR, UK) <giulio.caravagna@icr.ac.uk>\""

  list(REVOLVER_VERSION = REVOLVER_VERSION, REVOLVER_AUTHOR = REVOLVER_AUTHOR, REVOLVER_VERSION_NUMBER = REVOLVER_VERSION_NUMBER, REVOLVER_VERSION_CODNAME = REVOLVER_VERSION_CODNAME)
}




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
#
#
# obj_has_trees = function(x)
# {
#   if(is.null(x$phylogenies))
#     stop('Cannot proceed without having computed the trees -- are you calling the right function?')
# }


obj_has_jackknife = function(x)
{
  if(is.null(x$jackknife))
    stop('Cannot proceed without having computed the jackknife -- are you calling the right function?')
}


#' @importFrom utils object.size
saveFile = function(descr, fname, ...)
{
  input_list <- as.list(substitute(list(...)))
  save(input_list, file = fname)

  if(!is.na(descr))
  {
    pio::pioHdr("REVOLVER Save file", toPrint = c(`What`= descr), suffix = '\t -')

    stat = file.info(fname, extra_cols = FALSE)
    stat$size = format(utils::object.size(stat$size), "auto")
    stat$isdir = NULL
    stat$ctime = stat$atime = NULL

    print(stat)
    pio::pioStr("Path:", getwd(), suffix = '\n')
  }

}


########### Utilities to query a model



#### Functions to switch the internal representation of a model to an adjecaceny matrix, a data frame or other more compact notations


