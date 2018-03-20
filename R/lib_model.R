########### Utilities to query a model

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

