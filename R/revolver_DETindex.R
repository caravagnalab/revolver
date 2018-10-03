#' helper function to calculate diversity indices
#' 
#' @param   x           revolver output 
#' @param   n.boot      number of bootstraps (default is 1)
#' @param   table       return a table of results? (default is FALSE) 
#' @param   index       diversity index (Shannon or VA, default is Shannon) 
#' @param   type        'before.expansion' or 'after.expansion' (default after)
#' @param   driver.only use driver events only? (default is FALSE)
#' 
#' @return
#' 
#' @import  vegan
#' @import  qualvar
revolver_DETindex <- function(x, 
                              n.boot=1, 
                              table=FALSE, 
                              index='Shannon', 
                              type='after.expansion',
                              driver.only=FALSE) {
    
  feat <- revolver.featureMatrix(x)
  
  # switch result based on args 
  how <- paste0(ifelse(driver.only, "driver", "all"), 
                ifelse(type=="before.expansion", "Before", "After"))
  # message(how)
  res <- switch(how, 
    "driverBefore"=doidx(x, x$fit$multinomial.penalty, n.boot, table),
    "allBefore"=idx(x, x$fit$multinomial.penalty, n.boot, table),
    "driverAfter"=doidx(x, DFW2Matrix(feat$consensus.explosion), n.boot, table),
    "allAfter"=idx(x, DFW2Matrix(feat$consensus.explosion), n.boot, table)
  )
  return(res)


  # helper fn: Support 2 types of indices
  eval <- function(M, index) { 
    # {{{ Shannon or VA
    if (index == 'Shannon') {
      values <- vegan::diversity(M,MARGIN=2)/log(vegan::specnumber(M,MARGIN=2))
      values[is.na(values)] = 0 # normalized at 0 (100% homogenous trajectories)
    } else if (index == 'VA') {
      values <- apply(M, 2, qualvar::VA)
    } else { 
      stop('Unknown index? Use "Shannon" or "VA".')
    }
    return(values)
    # }}}
  } 

  # helper fn: all lesions index 
  idx <- function(x, w, n, table = F) { 
    # {{{ progress bar provided for indexing
    w <- w[, -1]
    DET <- eval(w, index) # empirical
    # DET = 0
    npbs <- c()
    pb <- txtProgressBar(min = 1, max = n.boot, style = 3)
    pb.status <- getOption('revolver.progressBar', default = TRUE)
    for (i in 1:n.boot) {
      # update progress bar
      if (pb.status) setTxtProgressBar(pb, i)
      npb <- w[, sample(ncol(w), replace = TRUE)]
      npbs[i] <- sum(eval(npb, index))
      # DET = DET + eval(npb, index)/n.boot
    }
    close(pb)
    cat(green(' DONE\n\n'))
    return(list(DET.cohort = npbs, DET.driver = DET))
    # }}}
  } 

  # helper fn: driver-only index 
  doidx <- function(x, w, n, table = F) {
    # {{{ driver-only empirical diversity
    w <- w[, -1]
    DET <- eval(w, index) # empirical
    res <- list(DET.cohort = NULL, DET.driver = DET)
    return(res) 
    # }}}
  }

}
