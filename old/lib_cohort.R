





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



