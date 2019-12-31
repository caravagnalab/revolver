# save loading of quartz for windows

.onLoad <- function(libname, pkgname) {
  
  options(pio.string_fg_colour = crayon::bgYellow$black)
  
  # =-=-=-=-=-=-
  # Required packages will be listed here
  # =-=-=-=-=-=-
  requirements = c('tidyverse', 'pio', 'crayon', 'ggpubr', 'RColorBrewer', 'cowplot',
                   'ctree', 'mtree')
  
  suppressMessages(sapply(requirements, require, character.only = TRUE))
  
  # =-=-=-=-=-=-
  # Package options
  # =-=-=-=-=-=-
  options(pio.string_fg_colour = crayon::bgYellow$black)
  
  # =-=-=-=-=-=-
  # Header
  # =-=-=-=-=-=-
  
  revolver_welcome_message =  getOption('revolver_welcome_message', default = TRUE)
  
  if(revolver_welcome_message)
  {
    pio::pioHdr('REVOLVER - Repeated Evolution in Cancer')
    pio::pioStr("Author : ", "Giulio Caravagna <gcaravagn@gmail.com>", suffix = '\n')
    pio::pioStr("GitHub : ", "caravagn/revolver", suffix = '\n')
    pio::pioStr("   WWW : ", "https://caravagn.github.io/revolver/", suffix = '\n')
    
    
    cat(
      "\n > REVOLVER is part of the", crayon::green("\"evoverse\""),
      crayon::blue("[https://bit.ly/2orn94e]"),
      "- a collection of packages to implement Cancer Evolution analyses from cancer sequencing data.\n"
    )
    
    
    options(revolver_welcome_message = FALSE)
  }
  
  invisible()
}

