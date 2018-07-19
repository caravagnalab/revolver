revolver_DETindex <-
  function(x,
           n.boot = 1,
           table = F,
           index = 'Shannon',
           type = 'after.expansion',
           driver.only = FALSE)
  {
    # require(qualvar)
    # require(vegan)

    # Support 2 types of indexes
    eval = function(M, index) {
      if (index == 'Shannon') {
        values = vegan::diversity(M, MARGIN = 2) / log(vegan::specnumber(M, MARGIN = 2))
        values[is.na(values)] = 0 # normalized when 0 (perfectly homogenous trajectories)
        return(values)
      }

      if (index == 'VA') {
        values = apply(M, 2, qualvar::VA)
        return(values)
      }

      stop('Unknown index? Use "Shannon" or "VA".')
    }

    indexFun = function(x, w, n, table = F)
    {
      w <- w[, -1]

      DET = eval(w, index) # empirical

      # DET = 0
      npbs <- c()

      pb = txtProgressBar(min = 1,
                          max = n.boot,
                          style = 3)
      pb.status = getOption('revolver.progressBar', default = TRUE)

      for (i in 1:n.boot) {
        # update progress bar
        if (pb.status)
          setTxtProgressBar(pb, i)

        npb = w[, sample(ncol(w), replace = TRUE)]

        npbs[i] <- sum(eval(npb, index))
        # DET = DET + eval(npb, index)/n.boot
      }

      close(pb)

      cat(green(' DONE\n\n'))
      return(list(DET.cohort = npbs, DET.driver = DET))
    }

    dr.indexFun = function(x, w, n, table = F)
    {
      w <- w[, -1]

      DET = eval(w, index) # empirical

      return(list(DET.cohort = NULL, DET.driver = DET))
    }


    feat = revolver.featureMatrix(x)

    if(driver.only) {

      if (type == 'before.expansion')
        return(dr.indexFun(x, w = x$fit$multinomial.penalty, n.boot, table))
      else
        return(dr.indexFun(x, w = DFW2Matrix(feat$consensus.explosion), n.boot, table))
    }


    if (type == 'before.expansion')
      return(indexFun(x, w = x$fit$multinomial.penalty, n.boot, table))
    else
      return(indexFun(x, w = DFW2Matrix(feat$consensus.explosion), n.boot, table))

    # results = list(
    #   `before` = indexFun(x, w = x$fit$multinomial.penalty, n.boot, table),
    #   `after` = indexFun(x, w = DFW2Matrix(feat$consensus.explosion), n.boot, table)
    # )
    #
    # return(results)
  }



#' Dump statistics for fit to an Excel file.
#'
#' @param x REVOLVER cohort fit
#' @param file output file
#'
#' @return none
#' @export
#' 
#' @importFrom pio pioHdr pioTit
#'
#' @examples
#' data(CRC.cohort)
#' fit = revolver_fit(CRC.cohort)
#' revolver_dumpStatistics(CRC.cohort)
revolver_dumpStatistics = function(x, file = 'REVOLVER-Statistics.xslx') {
  pioHdr('REVOLVER Dump statistics to Excel',
         "File: ", file,
         suffix = '\t')

  # require(xlsx)
  # dyn.load('/Library/Java/JavaVirtualMachines/jdk1.8.0_65.jdk/Contents/Home/jre/lib/server/libjvm.dylib')
  # require(rJava)

  pioTit('Features: drivers, occurrences & information transfer')

  ft = revolver.featureMatrix(x)

  xlsx::write.xlsx(
    x$dataset[x$dataset$is.driver,],
    file = file,
    sheetName = 'Drivers',
    row.names = FALSE
  )
  xlsx::write.xlsx(
    ft$occurrences,
    file = file,
    sheetName = 'Occurrences',
    append = TRUE
  )
  xlsx::write.xlsx(
    ft$consensus.explosion,
    file = file,
    sheetName = 'Information_Transfer',
    append = TRUE,
    row.names = FALSE
  )

  if (!is.null(x$fit)) {
    pioTit('Fit: multinomial, penalty & fit')

    xlsx::write.xlsx(
      x$fit$multinomial.penalty,
      file = file,
      sheetName = 'Multinomial',
      append = TRUE
    )
    xlsx::write.xlsx(
      x$fit$penalty,
      file = file,
      sheetName = 'Penalty',
      append = TRUE
    )
    xlsx::write.xlsx(
      x$fit$solutionID,
      file = file,
      sheetName = 'Fit',
      append = TRUE
    )
  }

  if (!is.null(x$cluster)) {
    pioTit('Clustering: distance, clusters & features per cluster')

    xlsx::write.xlsx(
      x$cluster$distances,
      file = file,
      sheetName = 'Distance',
      append = TRUE
    )
    xlsx::write.xlsx(
      x$cluster$clusters,
      file = file,
      sheetName = 'Clusters',
      append = TRUE
    )

    groups = sort(x$cluster$clusters)

    for (g in sort(names(x$cluster$labels.colors)))
    {
      ft = revolver.featureMatrix(x, patients = names(groups[groups == as.character(g)]))
      xlsx::write.xlsx(
        ft$consensus.explosion,
        file = file,
        sheetName = g,
        append = TRUE,
        row.names = FALSE
      )
    }
  }

  if (!is.null(x$jackknife)) {
    pioTit('Jackknife: co-clustering, clusters & features per cluster')

    xlsx::write.xlsx(
      x$jackknife$cluster,
      file = file,
      sheetName = 'Jackknife-clustering',
      append = TRUE
    )
    xlsx::write.xlsx(
      x$jackknife$edges,
      file = file,
      sheetName = 'Jackknife-edges-detection',
      append = TRUE,
      row.names = FALSE
    )
  }
}
