filter.clusters = function(d, cutoff.numMuts = 5)
{
  d = split(d, f = d$cluster)
  counts = lapply(d, nrow)

  # removed = d[counts <= cutoff.numMuts]
  # #
  # lapply(removed, function(w) print())

  return(
    list(
      removed = Reduce(rbind, d[counts <= cutoff.numMuts]),
      remained = Reduce(rbind, d[counts > cutoff.numMuts])
    ))

  # return(Reduce(rbind, d[counts > cutoff.numMuts]))
}

clusters.table = function(d, sample.groups)
{
  d = split(d, f = d$cluster)
  clusters = names(d)

  counts = lapply(d, nrow)

  CCF.means = lapply(d, function(x) {
    matrixStats::colMedians(as.matrix(x[, sample.groups, drop = F]))
  })
  CCF.means = Reduce(rbind, CCF.means)

  if(!is.matrix(CCF.means)) CCF.means = matrix(CCF.means, nrow = 1)
  colnames(CCF.means) = sample.groups

  R = data.frame(
    is.driver = unlist(lapply(d, function(x) any(x$is.driver))),
    is.clonal = unlist(lapply(d, function(x) any(x$is.clonal))),
    nMuts = unlist(counts),
    CCF.means
  )
  rownames(R) = clusters

  R = R[order(R$nMuts, decreasing = T), , drop = F]

  return(R)
}

nclusters = function(d){return(length(unique(d$cluster)))}

permuteClusterIds = function(d){
  ids = 1:length(unique(d$cluster.tracerx))
  names(ids) = unique(d$cluster.tracerx)
  d$cluster = ids[as.character(d$cluster.tracerx)]

  return(d[order(d$cluster), ])
}

panelToList = function(d)
{
  d = d[!d$excluded, ]
  d = d[d$parent != -1, ]
  d = d[, c('lab', 'parent')]
  colnames(d) = c('to', 'from')
  return(d)
}

hashTrees = function(clonevol.obj, sample.groups)
{
  trees = map.trees = NULL
  for(s in sample.groups)
  {
    l = lapply(
      clonevol.obj$models[[s]],
      function(x, region)
      {
        # hash ;)
        model = DataFrameToMatrix(x)
        model = sort(MatrixToEdges(model))

        return(list(model = paste(model, collapse = ":"), region = region))
      },
      region = s)

    trees = rbind(trees, Reduce(rbind, l))
  }

  rownames(trees) = NULL
  trees = data.frame(
    model = unlist(trees[, 'model']),
    model = unlist(trees[, 'region'])
  )
  colnames(trees) = c('model', 'region')

  # cat('Discarded empty trees :', length(which(trees$model == "")), '\n')
  trees = trees[trees$model != "", , drop = FALSE]

  # cat('Total trees  :', nrow(trees), '\n')
  # cat('Unique trees : ', length(unique(trees$model)), '\n')
  cat("Hashed", length(unique(trees$model)), 'trees\n')

  trees$model = as.factor(trees$model)

  map.trees = split(trees, f = trees$model)
  map.trees = lapply(map.trees, function(x)paste(x[, 'region'], collapse = ','))

  return(map.trees[map.trees != ""])

}


plotCCFTrees = function(map, my.data, binary.data, consensus, MI.table, distribution, distribution.scores, file, cex = 1, max.trees.in.plots = 50)
{
  ##### Plot the distribution
  ndistribution = ND = length(distribution)
  cat('* Plotting', ndistribution, 'trees for distribution\n')
  if(ndistribution > max.trees.in.plots){
    cat('\t Too many trees -- will plot only ', max.trees.in.plots, 'ranked by score..\n')

    ndistribution = max.trees.in.plots
    distribution = distribution[1:max.trees.in.plots]
  }

  toMerge = sapply(
    1:ndistribution,
    function(x){
      f = paste(file, x, 'distribution.pdf', sep = '-')

      plot(
        distribution[[x]],
        plot.stat = TRUE,
        graph.cex = cex,
        table.cex = cex,
        edge.label = round(MI.table, 2),
        file = f
        )

      return(f)
      })
  print(toMerge)

  rows = round(length(distribution)/5)
  cols = ceiling((length(distribution))/rows)
  cmd = paste('pdfjam ',
              paste(toMerge, collapse =' '),
              ' --nup ',  2, 'x', 2,  ' --landscape --outfile ', PATIENT, '-distribution.pdf', sep = '')
  print(cmd)
  system(cmd)
}

consensusModel = function(clonevol.obj, sample.groups)
{
  tr = unlist(lapply(clonevol.obj$models, length))
  # cat('* consensusModel over the following trees:', tr, ' -- TOT = ', sum(tr), '\n')

  S = lapply(clonevol.obj$models, function(x) Reduce(rbind, x))
  S = lapply(S, unique)
  S = Reduce(rbind, S)

  # S = Reduce(rbind, lapply(clonevol.obj$models, function(x) Reduce(rbind, x)))

  counts = DataFrameToEdges(S)
  counts = table(counts)

  counts = cbind(edgesToDataFrame(names(counts)), unlist(counts))
  # print(counts)
  # cat('* Consensus counts', counts$Freq, '\n')

  counts$counts = NULL

  counts = split(counts, f = counts$to)
  counts = lapply(counts, function(x) {y = x$Freq/ sum(x$Freq); names(y) = x$from; return(y)} )

  S = DataFrameToMatrix(S)
  S = MatrixToDataFrame(S)

  return(list(S = S, weights = counts))
}


all.possible.trees = function(
  G,
  W,
  sspace.cutoff = 10000,
  n.sampling = 1000
  )
{
  M = DataFrameToMatrix(G)
  r = root(M)
  parents = sapply(colnames(M), pi, model = M)

  parents = parents[!unlist(lapply(parents, is.null))]
  # cat("* Nodes: ", colnames(M), '\n')
  # cat("* Root: ", r, '\n')
  # cat("* Parent set \n\t")
  #
  # nothing = lapply(parents, function(x) cat(paste('{', paste(x, collapse =', ', sep =''), '}', collapse = '')))

  # singletons -- template
  singletons = parents[unlist(lapply(parents, function(x) length(x) == 1))]
  nsingletons = length(singletons)


  # cat('\n* Singletons:', names(singletons))
  # cat(' [ n =', nsingletons, ']\n')

  sglt = expand.grid(singletons, stringsAsFactors = FALSE)
  if(ncol(sglt) == 0) {
    sglt = NULL
  }

  alternatives = parents[unlist(lapply(parents, function(x) length(x) > 1))]
  nalternatives = length(alternatives)

  altn = NULL
  if(nalternatives > 0)
  {
    combalternatives = prod(unlist(lapply(alternatives, length)))
    # cat('* Alternatives:', names(alternatives), '-- num.', combalternatives, '\n')
    cat(combalternatives, 'solutions; ')


    if(combalternatives > sspace.cutoff)
    {
      cat("sampling is", red('Montecarlo for', n.sampling, 'distinct trees\n'))

      return(
        weighted.sampling(
          DataFrameToMatrix(G),
          W,
          n.sampling
        )
      )
    }
    else cat("sampling is", green('Exhaustive.\n'))




    # cat('* Generating and checking trees in 2 seconds ...\n')
    # Sys.sleep(2)

    altn = expand.grid(alternatives, stringsAsFactors = FALSE)
  }
  else cat(red('There are no alternatives!\n'))

  # all combinations
  if(is.null(altn) && is.null(sglt)) stop('Error -- no trees?')

  if(is.null(sglt) && !is.null(altn)) comb = altn
  if(!is.null(sglt) && is.null(altn)) comb = sglt
  if(!is.null(sglt) && !is.null(altn)) comb =  cbind(altn, sglt)


  pb = txtProgressBar(min = 1,
                      max = nrow(comb),
                      style = 3)
  pb.status = getOption('revolver.progressBar', default = TRUE)

  models = NULL
  for(i in 1:nrow(comb))
  {
    if (pb.status)
      setTxtProgressBar(pb, i)

    tree = data.frame(from = unlist(comb[i, ]), to = colnames(comb), stringsAsFactors = FALSE)
    test.tree = DataFrameToMatrix(tree)

    if(length(root(test.tree)) > 1 ) {
      # cat('Solution', DataFrameToEdges(tree), 'has multiple roots, removed\n')
      next;
    }

    if(!igraph::is_dag(igraph::graph_from_adjacency_matrix(test.tree))){
      # cat('Solution', DataFrameToEdges(tree), 'has loops, removed\n')
      next;
    }

    # # revert mapping
    # tree = DataFrameToMatrix(tree)
    # tree = reverse.mapping(tree, clusterIdsMapping)
    # tree = MatrixToDataFrame(tree)

    models = append(models, list(tree))
  }

  close(pb)


  # cat(length(models), ' ')
  # cat('\r')

  return(models)
}

binarize = function(d, sample.groups, cutoff = 1e-3)
{
  clusters = clusters.table(d, sample.groups)
  data = clusters[, sample.groups, drop = FALSE]
  data[data > cutoff] = 1
  data[data <= cutoff] = 0

  return(t(data))
}

reverse.mapping = function(M, map, from = rownames(M), cols = TRUE, rows = TRUE)
{
  new.rownames = NULL
  for(r in from)
    new.rownames = c(new.rownames, unique(map[map$cluster == r, 'cluster.tracerx']))

  if(cols) colnames(M) = new.rownames
  if(rows) rownames(M) = new.rownames

  return(M)
}


rankTrees = function(TREES, MI.table, structural.score)
{
  # cat('* Computing score for', length(TREES),'trees\n')



  pb = txtProgressBar(min = 1,
                      max = length(TREES),
                      style = 3)
  pb.status = getOption('revolver.progressBar', default = TRUE)

  MI.TREES = NULL
  for(i in 1:length(TREES))
  {
      if (pb.status) setTxtProgressBar(pb, i)

        M = DataFrameToMatrix(TREES[[i]])
        M = M[colnames(MI.table), colnames(MI.table)]

        M.entries = MI.table[which(M > 0, arr.ind = TRUE)]

        val = NA
        if(all(is.null(structural.score))) val = prod(M.entries)
        else val = prod(M.entries) * structural.score[i]

        MI.TREES = c(MI.TREES, val)
  }
  close(pb)

  # MI.TREES = sapply(
  #
  #   function(x)
  #   {
  #     # if(pb.status) cat('@ ', x, '\r')
  #     if (pb.status) setTxtProgressBar(pb, i)
  #
  #     M = DataFrameToMatrix(TREES[[x]])
  #     M = M[colnames(MI.table), colnames(MI.table)]
  #
  #     # print(TREES[[x]])
  #     # print(M)
  #     # print(structural.score[x])
  #     M.entries = MI.table[which(M > 0, arr.ind = TRUE)]
  #     # print(M.entries)
  #
  #     # M = hadamard.prod(M, MI.table)
  #     if(all(is.null(structural.score))) prod(M.entries)
  #     else prod(M.entries) * structural.score[x]
  #   })


  cat('Scores range from ', max(MI.TREES), '(max) to', min(MI.TREES), '(min)\n')
  # print(head(sort(table(MI.TREES), decreasing = T)))

  o = order(MI.TREES, decreasing = TRUE)
  MI.TREES = MI.TREES[o]
  TREES = TREES[o]
  structural.score = structural.score[o]

  return(list(TREES = TREES, SCORES = MI.TREES, PENALTIES = structural.score))
}


useClonevo = function(my.data, sample.groups, clonal.cluster)
{
  for(s in sample.groups)
  {
    if(max(my.data[, s]) != max(my.data[my.data$cluster == clonal.cluster, s])) {
      my.data[my.data$cluster == clonal.cluster, s] = 100
        # max(my.data[, s]) + 1
      warning('CCF correction for clonal cluster.')
    }
  }

  # Clonevo wants progressive IDs for clusters...
  my.data$cluster.tracerx = my.data$cluster
  my.data = permuteClusterIds(my.data)

  # cat('Creating progressive IDs for Clonevo\n')
  # print(clusters.table(my.data, sample.groups))

  clusterIdsMapping = my.data[, c('cluster', 'cluster.tracerx')]
  clusterIdsMapping = unique(clusterIdsMapping)
  rownames(clusterIdsMapping) = clusterIdsMapping$cluster
  clusterIdsMapping$cluster = NULL


  # Clonevol -- modified...
  capture.output({
    clonevol.obj = infer.clonal.models(
      variants = my.data,
      cluster.col.name = 'cluster',
      # vaf.col.names = vaf.col.names,
      ccf.col.names = sample.groups,
      sample.names = sample.groups,
      cancer.initiation.model = 'monoclonal',
      # cancer.initiation.model = 'polyclonal',
      # subclonal.test = 'bootstrap',
      subclonal.test = 'none',
      subclonal.test.model = 'non-parametric',
      # subclonal.test.model = 'beta-binomial',
      num.boots = 1000,
      founding.cluster = clonal.cluster,
      cluster.center = 'median',
      ignore.clusters = NULL,
      clone.colors = NULL,
      min.cluster.vaf = 0.01,
      # min probability that CCF(clone) is non-negative
      sum.p = 0.05,
      # alpha level in confidence interval estimate for CCF(clone)
      alpha = 0.05,
      verbose = F
    )
  })

  clonevol.obj$matched = NULL
  clonevol.obj$params = NULL
  clonevol.obj$variants = NULL
  clonevol.obj$num.matched.models = NULL

  for(s in sample.groups)
  {

    w = clonevol.obj$models[[s]]


    w = lapply(w, function(x){
      x = x[, c('lab', 'parent')]
      x = x[!is.na(x$parent), ]
      x = x[ x$parent != -1, ]
      colnames(x) = c('to', 'from')
      x = x[, c('from', 'to')]
      x$from = clusterIdsMapping[x$from, ]
      x$to = clusterIdsMapping[x$to, ]
      return(x)
    })

    clonevol.obj$models[[s]] = w

  }

  return(clonevol.obj)
}

computeMI.table = function(binary.data, MI.Bayesian.prior = 0, add.control = FALSE)
{
  if(add.control) binary.data = rbind(binary.data, wt = 0)

  # • a=0:maximum likelihood estimator (see entropy.empirical)
  # • a=1/2:Jeffreys’ prior; Krichevsky-Trovimov (1991) entropy estimator
  # • a=1:Laplace’s prior
  # • a=1/length(y):Schurmann-Grassberger (1996) entropy estimator
  # • a=sqrt(sum(y))/length(y):minimax prior

  MI.table = matrix(
    apply(
      expand.grid(colnames(binary.data), colnames(binary.data)),
      1,
      function(x) {
        # Counting process.. with a Bayesian prior
        i = x[1]
        j = x[2]
        jo11 = (binary.data[, i] %*% binary.data[, j])/nrow(binary.data)
        jo10 = (binary.data[, i] %*% (1-binary.data[, j]))/nrow(binary.data)
        jo01 = ((1-binary.data[, i]) %*% binary.data[, j])/nrow(binary.data)
        jo00 = 1 - (jo10 + jo01 + jo11)
        entropy::mi.Dirichlet(matrix(c(jo11, jo10, jo01, jo00), nrow = 2), a = MI.Bayesian.prior)
      }),
    byrow = TRUE, ncol = ncol(binary.data))
  colnames(MI.table) = rownames(MI.table) = colnames(binary.data)

  return(MI.table)
}


weightMI.byMultinomial = function(MI.table, W)
{
  Coeff = matrix(0, ncol = ncol(MI.table), nrow = nrow(MI.table))
  rownames(Coeff) = rownames(MI.table)
  colnames(Coeff) = colnames(MI.table)

  for(j in names(W))
    Coeff[ names(W[[j]]) , j] = unlist(W[[j]])

  return(matrixcalc::hadamard.prod(MI.table, Coeff))
}

weighted.sampling = function(G, W, n)
{
  S = S.hashcodes = NULL

  sampleT = function()
  {
    r = root(G)

    E = setdiff(colnames(G), r)
    E = sample(E, length(E))

    tree = NULL
    repeat {
      pi = sapply(E, function(node){
        parents = W[[node]]
        # print(node)
        # print(parents)
        draw = sample(names(parents), prob = unlist(parents), size = 1)
      })

      tree = data.frame(from = unlist(pi), to = names(pi), stringsAsFactors = FALSE)
      R = c(r, reach(tree, r))

      if(all(colnames(G) %in% R)) break
    }

    return(tree)
  }

  pb.status = getOption('revolver.progressBar', default = TRUE)

  c = 0
  repeat{
    Tree = sampleT()
    hash = paste(sort(DataFrameToEdges(Tree)), collapse = ':')

    if(!(hash %in% S.hashcodes))
    {
      S = append(S, list(Tree))
      S.hashcodes = c(S.hashcodes, hash)
      c = c + 1

      # cat('Found one tree')
    }
    # else {cat('Already sampled tree\n')}

    if(pb.status) cat('@ ', c, '\r')


    if(c == n) break;
  }

  #
  # cat('Sampled Trees: Cache (head)\n')
  # print(head(S.hashcodes))
  #
  return(S)

}

probability.raisers = function(G, W, n)
{
  S = S.hashcodes = NULL

  sampleT = function()
  {
    r = root(G)

    E = setdiff(colnames(G), r)
    E = sample(E, length(E))

    tree = NULL
    repeat {
      pi = sapply(E, function(node){
        parents = W[[node]]
        draw = sample(names(parents), prob = unlist(parents), size = 1)
      })

      tree = data.frame(from = unlist(pi), to = names(pi), stringsAsFactors = FALSE)
      R = c(r, reach(tree, r))

      if(all(colnames(G) %in% R)) break
    }

    return(tree)
  }

  c = 0
  repeat{
    Tree = sampleT()
    hash = paste(sort(DataFrameToEdges(Tree)), collapse = ':')

    if(!(hash %in% S.hashcodes))
    {
      S = append(S, list(Tree))
      S.hashcodes = c(S.hashcodes, hash)
      c = c + 1

      # cat('Found one tree')
    }
    # else {cat('Already sampled tree\n')}


    if(c == n) break;
  }


  cat('Sampled Trees: Cache (head)\n')
  print(head(S.hashcodes))

  return(S)
}


checkColinearity = function(TREES, ccf)
{
  muts = colnames(ccf)

  C = NULL
  for(m in muts)
  {
    conf = NULL
    for(n in muts){
      if((ccf[,m] %*% ccf[,n]) == 0) {
        cat('Conflict: ', m, ' and ',n, '\n')
        conf = c(conf, n)
      }
    }

    if(length(conf) > 0) {
      C = append(C, list(conf))
      names(C)[length(C)] = m
    }
  }

  cat('[check co-linearity] Violations\n')
  print(C)

  new.TREES = NULL
  for(t in 1:length(TREES))
  {
    TR = TREES[[t]]

    violates = FALSE
    for(m in muts){
      f = reach(TR, m)
      if(!is.null(f) && any(f %in% C[[m]]))
      {
        cat('Tree ', DataFrameToEdges(TR), 'violates co-linearity since', m, '-->*', intersect(f, C[[m]]), '\n')
        violates = TRUE

        # g = graph_from_adjacency_matrix(DataFrameToMatrix(TR))
        # plot(g)
        # stop('ss')
      }
    }

    if(!violates) new.TREES = append(new.TREES, list(TR))

  }

  return(new.TREES)
}


# for every edge  x --> y, the number of times that the CCF of x is greater than the CCF of y
edge.penalty.for.direction = function(TREES, ccf)
{
  nodes = rownames(ccf)

  p = expand.grid(nodes, nodes)
  colnames(p) = c('from', 'to')

  direction = apply(p, 1, function(x){
    1 - sum(as.numeric(ccf[x[1], ] < ccf[x[2], ]))/ncol(ccf)
  })
  p = cbind(p, direction)

  rownames(p) = DataFrameToEdges(p)

  cat('* Penalty table\n')
  print(p)

  cat('* Computing penalty for ', length(TREES),' trees\n.')

  scores = sapply(1:length(TREES), function(x)
  {
    cat('@ ', x, '\r')
    prod(p[DataFrameToEdges(TREES[[x]]), 'direction'])
  })

  return(scores)
}

# for every node  x --> y1 ... yK, the number of times that the CCF of x is greater than the sum of the CCFs of y1 ... yK
node.penalty.for.branching = function(TREES, ccf)
{
  nodes = rownames(ccf)

  if(length(TREES) > 1)
  {
    pb = txtProgressBar(min = 1,
                      max = length(TREES),
                      style = 3)
    pb.status = getOption('revolver.progressBar', default = TRUE)
  }

  scores = NULL
  for(x in 1:length(TREES))
  {
    # update progress bar
    if(length(TREES) > 1)
      if (pb.status) setTxtProgressBar(pb, x)

    t = DataFrameToMatrix(TREES[[x]])


    c = sapply(nodes, function(n) {
      cl = children(t, n)
      if(length(cl) == 0) return(1)


      1 - sum(as.numeric(ccf[n, ] < colSums(ccf[cl, , drop = FALSE])))/ncol(ccf)
    })

    scores = c(scores, prod(c))
  }

  if(length(TREES) > 1)
    close(pb)


  return(scores)
}


# Compute Suppes' poset
poset = function(x, regions)
{
  samples = x[, regions, drop = FALSE]
  variables = rownames(x)

  marginals = rowSums(samples)/length(regions)
  bmarginals = lapply(names(marginals), function(w) marginals >= marginals[as.character(w)])
  bmarginals = lapply(bmarginals, function(w) w[which(w)])
  names(bmarginals) = names(marginals)

  joints = outer(
    1:nrow(samples), 1:nrow(samples),
    FUN = Vectorize( function(i,j) sum(samples[i, ] * samples[j,])/length(regions) ))
  colnames(joints) = rownames(joints) = rownames(samples)

  bmarginals = lapply(names(bmarginals), function(b){
    to = b
    bm = bmarginals[[b]]
    bool = sapply(names(bm), function(from)
      {
        joints[as.character(to), as.character(from)] >= marginals[as.character(to)] * marginals[as.character(from)]
      })
    names(bool) = names(bm)
    bool
  })
  names(bmarginals) = names(marginals)
  for(b in names(bmarginals)) bmarginals[[b]][b] = FALSE

  bmarginals = lapply(bmarginals, function(w) w[which(w)])
  bmarginals = bmarginals[sapply(bmarginals, function(w) length(w) > 0)]


  return(bmarginals)
}


compute.clusters = function(cohort, PARSERFUN){
  cohort$cluster = NA
  cohort$CCF = as.character(cohort$CCF)

  bin2int <- function(x)
  {
    x <- as.character(as.numeric(x))
    b <- as.numeric(unlist(strsplit(x, "")))
    pow <- 2 ^ ((length(b) - 1):0)
    sum(pow[b == 1])
  }

  for(i in 1:nrow(cohort))
  {
    binary.region = as.numeric(PARSERFUN(cohort$CCF[i]))
    cohort$cluster[i] = bin2int(binary.region)
  }

  patients = unique(cohort$patientID)

  for (patient in patients)
  {
    sub.cohort = cohort[which(cohort$patientID==patient),]
    clust.patient = data.frame(match = unique(sub.cohort$cluster))
    cohort$cluster[which(cohort$patientID==patient)] = match(sub.cohort$cluster, clust.patient$match)
  }
  return (cohort)

}
