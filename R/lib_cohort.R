
rev_count_information_transfer_comb = function(x, p) {
  stopifnot(!is.null(x$phylogenies))
  stopifnot(p %in% names(x$phylogenies))

  keys = lapply(
    x$phylogenies[[p]],
        function(w)
          paste(sort(DataFrameToEdges(w$transfer$driver)), collapse = ' '))

  keys = Reduce(rbind, keys)
  length(unique(keys))
}




CCF.parser = function(x)
{
  tk = strsplit(x, ';')[[1]]
  tk = unlist(strsplit(tk, ':'))

  samples = tk[seq(1, length(tk), 2)]

  values = tk[seq(2, length(tk), 2)]
  names(values) = samples

  return(values)
}



clonal.subclonal.table = function(x)
{
  clonal.drivers = table(x$dataset[x$dataset$is.driver & x$dataset$is.clonal, 'variantID'])
  subclonal.drivers = table(x$dataset[x$dataset$is.driver & !x$dataset$is.clonal, 'variantID'])

  tableD = data.frame(
    variantID = unique(x$dataset[x$dataset$is.driver, 'variantID']),
    stringsAsFactors = FALSE)
  rownames(tableD) = tableD$variantID

  tableD = cbind(tableD, Clonal = 0)
  tableD = cbind(tableD, SubClonal = 0)

  for(d in tableD$variantID){
    tableD[d, 'Clonal'] = clonal.drivers[d]
    tableD[d, 'SubClonal'] = subclonal.drivers[d]
  }
  tableD[is.na(tableD)] = 0
  tableD$Counts = tableD$Clonal + tableD$SubClonal

  tableD = tableD[order(tableD$Counts, decreasing = T), ]
  tableD$variantID = NULL

  tableD
}


