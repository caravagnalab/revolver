# =-=-=-=-=-=-=-=-=-=-=-=-=-
# These function compute the information transfer from REVOLVER trees
# =-=-=-=-=-=-=-=-=-=-=-=-=-

# This function takes the list of drivers in x and traverses backward the tree
# to determine the transitive closure used by REVOLVER's algorithm
combination_of_information_transfer = function(x, patient)
{
  if(!has_patient_trees(x, patient)) return(0)

  keys = lapply(
    seq_along(x$phylogenies[[patient]]),
      function(w)
        paste(sort(ctree:::DataFrameToEdges(ITransfer(x, patient, type = 'clones', rank = w))), collapse = ' ')
    )

    keys = Reduce(rbind, keys)
    length(unique(keys))
}
