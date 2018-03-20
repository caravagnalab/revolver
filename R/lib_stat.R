revolver_DETindex <- function(x, n.boot = 1, table = F, index = 'Shannon')
{
  if(is.null(x$fit)) stop("Fit a model first?")

  # require(qualvar)
  # require(vegan)

  # Support 2 types of indexes
  eval = function(M, index) {
    if(index == 'Shannon'){
      values = vegan::diversity(M, MARGIN = 2)/log(vegan::specnumber(M, MARGIN = 2))
      values[is.na(values)] = 0 # normalized when 0 (perfectly homogenous trajectories)
      return(values)
    }

    if(index == 'VA'){
      values = apply(M, 2, VA::VA)
      return(values)
    }

    stop('Unknown index? Use "Shannon" or "VA".')
  }

  indexFun = function(x, w, n, table = F)
  {
    w <- w[,-1]

    DET = eval(w, index)
    npbs <- c()

    pb = txtProgressBar(min = 1, max = n.boot, style = 3)

    cat(cyan('Bootstrapping the DET index \n\tNon-parametric replicates n = ', n.boot,'\n\tIndex type: ', index, '\n'))

    for( i in 1:n.boot) {
      setTxtProgressBar(pb, i) # update progress bar

      npb = w[, sample(ncol(w), replace = TRUE)]
      npbs[i] <- sum(eval(npb, index))
    }

    cat(green(' DONE\n'))
    return(list(DET.cohort = npbs, DET.driver = DET))
  }

  feat = revolver.featureMatrix(x)
  results = list(
    `before` = indexFun(x, w = x$fit$multinomial.penalty, n.boot, table),
    `after` = indexFun(x, w = DFW2Matrix(feat$consensus.explosion), n.boot, table)
  )

  return(results)
}


revolver_dump_statistics = function(x, file = 'REVOLVER-Statistics.xslx') {
  # require(xlsx)
  # dyn.load('/Library/Java/JavaVirtualMachines/jdk1.8.0_65.jdk/Contents/Home/jre/lib/server/libjvm.dylib')
  # require(rJava)


  ft = revolver.featureMatrix(x)

  write.xlsx(x$dataset[fit$dataset$is.driver, ], file = file, sheetName = 'Drivers', row.names = FALSE)
  write.xlsx(ft$occurrences, file = file, sheetName = 'Occurrences', append = TRUE)
  write.xlsx(ft$consensus.explosion, file = file, sheetName = 'Information_Transfer', append = TRUE,  row.names = FALSE)

  if(!is.null(x$fit)) {
    write.xlsx(x$fit$multinomial.penalty, file = file, sheetName = 'Multinomial', append = TRUE)
    write.xlsx(x$fit$penalty, file = file, sheetName = 'Penalty', append = TRUE)
    write.xlsx(x$fit$solutionID, file = file, sheetName = 'Fit', append = TRUE)
  }

  if(!is.null(x$fit)) {
    write.xlsx(x$cluster$distances, file = file, sheetName = 'Distance', append = TRUE)
    write.xlsx(x$cluster$clusters, file = file, sheetName = 'Clusters', append = TRUE)

    groups = sort(x$cluster$clusters)

    for(g in sort(names(x$cluster$labels.colors)))
    {
      ft = revolver.featureMatrix(x, patients = names(groups[groups == as.character(g)]))
      write.xlsx(ft$consensus.explosion, file = file, sheetName = g, append = TRUE,  row.names = FALSE)
    }
  }

}
