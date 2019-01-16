CCF_phylogeny_univariate = function(x, clonal.cluster)
{
  # Compute all posssible branches out of s from x
  branching = function(
    n, # node
    x  # CCF values
  )
  {
    s = which(names(x) == n)
    y = x[-s]

    size = 1
    children = NULL
    repeat {
      ch = combn(names(y), size)
      sums = apply(ch, 2, function(w) sum(y[w]))

      ch = ch[, sums <= x[n], drop = FALSE]
      ch = apply(ch, 2, list)

      # print(ch)

      if(length(ch) == 0) break;

      children = append(children, ch)
      size = size + 1

      if(size > length(y)) break;
    }

    children
  }

  # Empty adjMat for vector of names x, with a fake WL node
  empty = function(x) {

    m = matrix(0, ncol = length(x), nrow = length(x))
    colnames(m) = rownames(m) = names(x)

    # WL
    m = cbind(m, `WL` = 0)
    m = rbind(m, `WL` = 0)

    # r = names(x)[which.max(x)]
    r = clonal.cluster
    m['WL', r] = 1

    m
  }

  # Queue functions
  pop = function(Q){ Q[-1] }
  push = function(Q,x){ append(Q, list(x)) }
  peek = function(Q){ Q[[1]] }

  # Queue with cache
  # push.cache = function(Q,x) {
  # sgn = lapply(Q, as.vector)
  # sgn = sapply(sgn, paste, collapse = '')
  #
  # key = paste(as.vector(x), collapse = '')
  #
  # if(!(key %in% sgn)) Q = append(Q, list(x))
  #
  # Q
  # }

  # cache for trees
  cacheify = function(Q, x) {
    key = paste(as.vector(x$tree), collapse = '')
    Q = c(Q, key)
    Q
  }

  not_in_cache = function(Q, x) {
    if(length(Q) == 0) return(TRUE)
    key = paste(as.vector(x$tree), collapse = '')
    !(key %in% Q)
  }

  # push.cache2 = function(Q,x) {
  #   sgn = lapply(Q, function(w) as.vector(w$tree))
  #   sgn = sapply(sgn, paste, collapse = '')
  #
  #   key = paste(as.vector(x$tree), collapse = '')
  #
  #   if(!(key %in% sgn)) Q = append(Q, list(x))
  #
  #   Q
  # }

  # struct to start with is a tree with
  # r = names(x)[which.max(x)]
  r = clonal.cluster
  tree = empty(x)

  init = list(
    tree = tree,
    leaves = r,
    missing = setdiff(names(x), r)
  )

  # Q is a the queue of partial trees to process
  # TR are the final set of trees to return
  # CACHE is the set of substractures that we have already included in Q (those avoid duplicates)
  Q = TR = CACHE_Q = CACHE_TR = NULL

  Q = push(Q, init)
  CACHE_Q = cacheify(CACHE_Q, init)

  repeat{
    if(length(Q) == 0) break;

    # cat("\n\nQ ", length(Q), " cache ",length(CACHE_Q))
    # cat("\nTR ", length(TR), " cache ",length(CACHE_TR), '\n')

    model = peek(Q)
    Q = pop(Q)

    # for every leaf node, get all its possible branches,
    # create corresponding trees and queue them
    for(n in model$leaves)
    {
      # take the list of nodes that are still missing,
      # and use it to compute the branches
      missing = x[c(n, model$missing)]
      branches = branching(n, missing)

      # for every branch,
      for(br in 1:length(branches))
      {
        br = unlist(branches[[br]])

        br.model = model
        br.model$leaves = c(setdiff(br.model$leaves, n), br)
        br.model$tree[n, br] = 1
        br.model$missing = setdiff(br.model$missing, br)

        # if there are no missing nodes, this is a final tree, otherwise it is
        # a partial and then it goes into Q. In both cases we avoid to include
        # somehting that we have already analyzed.
        if(length(br.model$missing) == 0)
        {
          if(not_in_cache(CACHE_TR, br.model))
          {
            TR = push(TR, br.model$tree)
            CACHE_TR = cacheify(CACHE_TR, br.model)
          }
          # else cat("- TR RJ")

        }
        else
        {
          if(not_in_cache(CACHE_Q, br.model))
          {
            Q = push(Q, br.model)
            CACHE_Q = cacheify(CACHE_Q, br.model)
          }
          # else cat("- Q RJ")
        }
      }
    }
  }

  TR
}

ClonEvol_surrogate = function(clusters, samples, clonal.cluster, min.CCF = 0.01)
{
  stopifnot(all(samples %in% colnames(clusters)))
  clusters = clusters[, samples, drop = FALSE]


  trees = lapply(samples,
                 function(w) {
                   regname = w

                   w = clusters[, w]
                   names(w) = rownames(clusters)

                   # CCF for this sample, we consider only clusters values above 1% CCF
                   w = w[w > min.CCF]

                   cat("\n Region: ", regname, " - #CCF clusters > 1%: ", length(w))

                   # tree -- default case with no tree
                   TR = list(data.frame(from = NULL, to = NULL))

                   # Check for consistency of the following parameters
                   # - the clonal cluster is always the maximum CCF value in all biopsies
                   # - there are no samples with all CCF below
                   if(length(w) > 1)
                   {
                     # check the clonal cluster
                     subclonal.clusters = setdiff(names(w), clonal.cluster)

                     if(length(subclonal.clusters) > 0 & !all(w[clonal.cluster] >= w[subclonal.clusters])) {
                       warning(
                         '[CORRECTION] The clonal cluster has CCF ', w[clonal.cluster],
                         ', lower than a subclonal cluster; we set it to the max in the sample...')

                       w[clonal.cluster] = max(w) + 0.01
                     }

                     # trees from this sample, it is a list
                     TR = CCF_phylogeny_univariate(w, clonal.cluster)

                     # turn them into the format that the tools expects and remove the temporary WL node
                     TR = lapply(TR, revolver:::MatrixToDataFrame)
                     TR = lapply(TR, function(y){ y = y[y$from != 'WL', , drop = FALSE] })
                   }
                   else {
                     warning('[SKIP] One region does not have at least 2 CCF clusters above 1%, and will be skipped.')
                   }

                   TR
                 }
  )
  names(trees) = samples


  # just throw an error if no regions have  remove samples with no trees associated, and throw an error if there are no trees at all
  non_empty_trees = sapply(trees, function(w) { any(lapply(w, nrow) > 0) } )

  if(all(!non_empty_trees)) stop('This patient has no trees, raising an error. Check you CCF estimates ...')

  trees
}
