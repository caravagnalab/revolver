# a-la-dictionary
driver = function(x, c){ x$dataset[which(x$dataset$cluster == c & x$dataset$is.driver), 'variantID'] }
driver.indexOf = function(x, d){ x$dataset[which(x$dataset$variantID == d & x$dataset$is.driver), 'cluster'] }

# Compute the frontier of "var" in a model. The frontier is the set of
# mutations in Gamma that are directly reachable via
# the transitive closure of the --> relation. So, these are the events
# selected by "var"'s evolutionary trajectories.
information.transfer = function(x, transitive.closure = FALSE, indistinguishable = FALSE)
{
  aux = function(r)
  {
    r.d = driver(x, r)

    if(length(r.d) > 0) return(r) # stop recursion

    c = children(model, r)

    if(is.null(c)) return(NULL) # leaf

    # recursion, reduction etc.
    return(
      Reduce(union,
             sapply(c, aux))
    )
  }

  model = x$adj_mat
  nodes.drivers = as.character(unique(x$dataset[x$dataset$is.driver, 'cluster']))

  # Root
  df = expand.grid(from = x$root, to = aux(x$root), stringsAsFactors = FALSE)

  # and all other stuff
  for(n in nodes.drivers)
  {
    df.n = NULL

    for(c in children(model, n))
      df.n = rbind(df.n, expand.grid(from = n, to = aux(c), stringsAsFactors = FALSE))

    df = rbind(df, df.n)
  }

  # Then, for all clones, we expand the drivers that they have
  expanded = apply(
    df,
    1,
    function(w) {
      e.from = driver(x, w['from'])
      if(length(e.from) == 0) e.from = w['from']

      e.to = driver(x, w['to'])

      expand.grid(from = e.from, to = e.to, stringsAsFactors = FALSE)
    })
  expanded = Reduce(rbind, expanded)

  # Then we see if we want to expand also the indistinguishible
  if(indistinguishable){
    for(n in nodes.drivers){
      expanded = rbind(expanded,
        expand.grid(from = driver(x, n), to = driver(x, n), stringsAsFactors = FALSE))
    }
  }

  # If we need transitive closures, we compute them here
  if(transitive.closure)
  {
    df = DataFrameToMatrix(df)
    df = nem::transitive.closure(df, mat = T, loops = FALSE)

    expanded = DataFrameToMatrix(expanded)
    expanded = nem::transitive.closure(expanded, mat = T, loops = FALSE)

    df = MatrixToDataFrame(df)
    expanded = MatrixToDataFrame(expanded)
  }

  return(list(clones = df, drivers = expanded))
}


find.path = function(M, from, to)
{
  # ASSUMPTION: the path exists

  # reverse edges' direction
  M = t(M)
  S = NULL

  # swap from and to
  bfrom = to
  to = from
  from = bfrom

  repeat
  {
    n = children(M, from)

    S = rbind(S, c(from = from, to = n))
    from = n
    if(n == to) break
  }

  S = data.frame(from = S[, 2], to = S[, 1])
  return(S)
}


DFW2Matrix = function(df){
  entries.names = unique(unlist(df[ c(1,2)]))

  M = matrix(0, ncol = length(entries.names), nrow = length(entries.names))
  colnames(M) = rownames(M) = entries.names

  for(i in 1:nrow(df))
    M[df[i, 'from'], df[i, 'to']] = df[i, 'count']
  M
}


# Compute the root(s) of a model
root = function(model){
  s = colSums(model)
  return(names(s[s==0]))
}

# Compute the leaves of a model
leaves = function(model){
  s = rowSums(model)
  model = t(model)
  return(names(s[s==0]))
}

# Compute the children of "var" in a model
children = function(model, var)
{
  model = model[var, ]
  model = model[model == 1]

  if(is.null(model)) return(NULL)
  return(names(model))
}

# Compute the frontier of "var" in a model. The frontier is the set of
# mutations in Gamma that are directly reachable via
# the transitive closure of the --> relation. So, these are the events
# selected by "var"'s evolutionary trajectories.
frontier = function(model, var, SHARED.VARIABLES, DICTIONARY)
{
  aux = function(r)
  {
    r.d = dictionary.get(DICTIONARY, r)

    if(any(r.d %in% SHARED.VARIABLES)) return(r.d) # stop recursion

    c = children(model, r)

    if(is.null(c)) return(NULL) # leaf

    # recursion, reduction etc.
    return(
      Reduce(union,
             sapply(c, aux))
    )
  }

  return(aux(var))
}

# Return the parent of "variable" in model
pi = function(model, variable)
{
  model = model[, variable]
  if(any(model > 0)) model = model[model > 0]
  else return(NULL)
  return(names(model))
}


# Given a Data Frame representation of a model, compute the set of
# nodes reachable from x
reach = function(df, x)
{
  if(!any(df$from == x)) return(NULL)

  dfB = df[ df$from == x, , drop = F]
  r = dfB$to

  # print(dfB)
  for(i in 1:length(r))
  {
    r = c(r, reach(df, r[i]))
  }

  return(r)
}

