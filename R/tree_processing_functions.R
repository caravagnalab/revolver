# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Functions that are used to manipulate trees within REVOLVER
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Three possible representaitons for a tree:
# - as edges strings "A~B", "B~C" etc.
# - as dataframe with from/ to clumns
# - as an adjacency matrix
# We have marshalling functions to switch among these representations
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# A~B to a DataFrame (from/to)
edgesToDataFrame = function(edges)
{
  # We convert the list of edges to bnlearn's from/to format
  dfedges = data.frame(from = NULL, to = NULL, stringsAsFactors = F)

  if(length(edges) == 0) return(dfedges)

  for(j in 1:length(edges))
  {
    aux = strsplit(edges[j], '~')[[1]]
    dfedges = rbind(dfedges, data.frame(from = aux[1], to = aux[2], stringsAsFactors = FALSE))
  }
  return(dfedges)
}

# A~B to a Adjacency Matrix
edgesToMatrix = function(edges)
{
  df = edgesToDataFrame(edges)
  vars = unique(unlist(df))
  matrix = matrix(0, nrow = length(vars), ncol = length(vars))
  colnames(matrix) = vars
  rownames(matrix) = vars

  if(nrow(df) == 0) return(matrix)

  for(j in 1:nrow(df))
  {
    matrix[df[j, 'from'], df[j, 'to']] = 1
  }

  return(matrix)
}

# Adjacency Matrix to a DataFrame
MatrixToDataFrame = function(matr)
{
  dfedges = data.frame(stringsAsFactors = F)
  for(i in 1:nrow(matr)) {
    for(j in 1:ncol(matr)){
      if(matr[i,j] == 1)
        dfedges = rbind(dfedges, data.frame(
          from = rownames(matr)[i],
          to = colnames(matr)[j],
          stringsAsFactors = FALSE))

    }
  }
  return(dfedges)
}

# Adjacency Matrix to A~B
MatrixToEdges = function(matr){
  return(DataFrameToEdges(MatrixToDataFrame(matr)))
}

# DataFrame to A~B
DataFrameToEdges = function(edges)
{
  edg = NULL
  for(j in 1:nrow(edges))
    edg = c(edg, paste(edges[j, 'from'], edges[j, 'to'], sep = '~'))
  return(edg)
}

# DataFrame to Adjacency Matrix
DataFrameToMatrix = function(edges)
{
  return(edgesToMatrix(DataFrameToEdges(edges)))
}

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Functions to query the properties of a tree, to check its
# consistency etc. All these functions assume to be working
# with the adjacency matrix representation of a tree
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# Return true if M is a tree
is_tree = function(M)
{
  cS = colSums(M)

  # - No more than 1 parent
  # - 1 root
  # - No nodes with 0 parents
  (max(cS) == 1) && (length(root(M)) == 1) && (sum(cS) == length(cS) - 1)
}

# Compute the root(s) of a model
root = function(model){
  s = colSums(model)
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


