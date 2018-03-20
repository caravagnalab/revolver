# Setup parallel for multi-core executions
setup_parallel = function(cores.ratio)
{
  # require(parallel)
  # require(doParallel)
  # require(crayon)

  # set the number of cores to be used in the parallelization
  cores = as.integer(cores.ratio * (detectCores() - 1))
  if (cores < 1) cores = 1

  cat(cyan('Registering to use', cores, 'cores out of', detectCores(), 'via \"parallel\" ...'))
  cl = makeCluster(cores)
  registerDoParallel(cl)
  cat(bgGreen(" OK.\n"))

  return(cl)
}

# Stop parallel for multi-core executions
stop_parallel =  function(cl)
{
  cat(cyan('Stopping parallel clusters ...'))
  stopCluster(cl)
  cat(bgGreen(" OK.\n"))
}

