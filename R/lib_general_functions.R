# nmfy = function(...)
# {
#   l = list(...)
#   # print(l)
#   print(match.call())
#   # print(match.arg(arg = '...', choices = '...'))
# }

revolver_version = function()
{
  REVOLVER_VERSION_NUMBER = "1.0"
  REVOLVER_VERSION_CODNAME = "\"Haggis and tatties\""

  REVOLVER_VERSION = paste(REVOLVER_VERSION_NUMBER, REVOLVER_VERSION_CODNAME, sep = ' ~ ')

  REVOLVER_AUTHOR = "\"Giulio Caravagna (ICR, UK) <giulio.caravagna@icr.ac.uk>\""

  list(REVOLVER_VERSION = REVOLVER_VERSION, REVOLVER_AUTHOR = REVOLVER_AUTHOR, REVOLVER_VERSION_NUMBER = REVOLVER_VERSION_NUMBER, REVOLVER_VERSION_CODNAME = REVOLVER_VERSION_CODNAME)
}


obj_has_clusters = function(x)
{
  if(is.null(x$cluster))
    stop('Cannot proceed without having computed clusters -- are you calling the right function?')
}

obj_has_fit = function(x)
{
  if(is.null(x$fit))
    stop('Cannot proceed without having computed fit models -- are you calling the right function?')
}

obj_has_evodistance = function(x)
{
  if(is.null(x$cluster$distances))
    stop('Cannot proceed without having computed the evolutionary distance -- are you calling the right function?')
}


obj_has_trees = function(x)
{
  if(is.null(x$phylogenies))
    stop('Cannot proceed without having computed the trees -- are you calling the right function?')
}


obj_has_jackknife = function(x)
{
  if(is.null(x$jackknife))
    stop('Cannot proceed without having computed the jackknife -- are you calling the right function?')
}


saveFile = function(descr, fname, ...)
{
  save(..., file = fname)

  if(!is.na(descr))
  {
    cat(crayon::cyan(descr), '\n')

    stat = file.info(fname, extra_cols = FALSE)
    stat$size = utils:::format.object_size(stat$size, "auto")
    stat$isdir = NULL
    stat$ctime = stat$atime = NULL

    print(stat)
    cat(crayon::cyan("Path:"), getwd(), '\n')
  }

}


########### Utilities to query a model



#### Functions to switch the internal representation of a model to an adjecaceny matrix, a data frame or other more compact notations

# A~B to a DataFrame (from/to)
edgesToDataFrame = function(edges)
{
  # We convert the list of edges to bnlearn's from/to format
  dfedges = data.frame(from = NULL, to = NULL, stringsAsFactors = F)
  # colnames(dfedges) = c('from', 'to')

  # print(dfedges)
  # print(length(edges))
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
