# save loading of quartz for windows

.onLoad <- function(libname, pkgname) {
  
  options(pio.string_fg_colour = crayon::bgYellow$black)
  
  suppressMessages(library(tidyverse))
  suppressMessages(library(dplyr))
  
  pio::pioHdr('REVOLVER - Repeated Evolution in Cancer')
  pio::pioStr("Author : ", "Giulio Caravagna <gcaravagn@gmail.com>", suffix = '\n')
  pio::pioStr("GitHub : ", "caravagn/revolver", suffix = '\n')
  
  
  # No longer need
    # if(.Platform$OS.type=="windows") {
    #     quartz<-function(...) windows(...)
    # }
}

