# save loading of quartz for windows

.onLoad <- function(libname, pkgname) {
    if(.Platform$OS.type=="windows") {
        quartz<-function() windows()
    }
}

