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
