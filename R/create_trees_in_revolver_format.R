create_trees_in_revolver_format = function(options, TREES, SCORES, patient, x, samples)
{
  cat("Total number of trees with non-zero scores is", length(TREES))
  
  if(length(TREES) > options['store.max'])
  {
    cat(red(' -- too many, storing top', options['store.max'], '\n'))
    
    TREES = TREES[1:as.numeric(options['store.max'])]
    SCORES = SCORES[1:as.numeric(options['store.max'])]
    
  } else cat(green(' -- storing all\n'))
  
  LSCORES = as.data.frame(SCORES)
  LSCORES = split(LSCORES, f = LSCORES[,1])
  
  LSCORES = lapply(LSCORES, function(w) w[sample(1:nrow(w), replace = FALSE), , drop = FALSE])
  
  # Shuffle indexes of equal-scoring fits to avoid a sistematic bias in the way I coded this
  permuted.indexes = as.integer(rev(unlist(lapply(LSCORES, rownames))))
  names(permuted.indexes) = NULL
  
  TREES = TREES[permuted.indexes]
  SCORES = SCORES[permuted.indexes]
  
  #####################################################################
  pb = txtProgressBar(min = 0,
                      max = length(TREES),
                      style = 3)
  pb.status = getOption('revolver.progressBar', default = TRUE)
  
  REVOLVER.TREES = NULL
  for(i in 1:length(TREES))
  {
    if (pb.status)
      setTxtProgressBar(pb, i)
    
    tree = revolver_phylogeny(
      x,
      patient = patient,
      M = TREES[[i]],
      score = SCORES[i],
      annotation = paste('Ranked ', i, '/', length(TREES), sep ='')
    )
    
    REVOLVER.TREES = append(REVOLVER.TREES, list(tree))
  }
  
  close(pb)
  #####################################################################
  
  
  
  REVOLVER.TREES
}
