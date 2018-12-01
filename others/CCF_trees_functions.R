#
# # Function that replaces all we need from clonEvol. It creates
# # for a single sample (CCF values) all possible trees that satisfy the
# # pigeonhole principle
# single_sample_CCF_pigeonhole_generator = function(m, cutoff)
# {
#   all.trees = NULL
#
#   for(i in 1:ncol(m))
#   {
#     pio::pioStr("\nAnalyzing sample CCF for", colnames(m)[i], suffix = '\n')
#
#     df = m[, i]
#     names(df) = rownames(m)
#     df = df[df >= cutoff]
#
#     edges = lapply(names(df), function(w) data.frame(from = w, to = names(df[df[w] >= df])))
#     edges = revolver:::DataFrameToMatrix(Reduce(rbind, edges))
#     diag(edges) = 0
#
#     edges = revolver:::MatrixToDataFrame(edges)
#     # plot(graph_from_edgelist(as.matrix(edges)))
#
#     span_graphs = combn(1:nrow(edges), length(df) - 1)
#     span_trees = NULL
#
#     pb = txtProgressBar(0, ncol(span_graphs), style = 3)
#
#     for(j in 1:ncol(span_graphs))
#     {
#       setTxtProgressBar(pb, j)
#
#       tree = edges[span_graphs[, j], ]
#       tree = revolver:::DataFrameToMatrix(tree)
#
#       # Structural constraint, it should have no
#       # - multiple roots (disconnected components)
#       # - no root
#       # - loops
#       # - a node with multiple parents
#
#       multi_root = sum(apply(tree, 2, function(w) all(w == 0))) > 1
#       no_root = sum(apply(tree, 2, function(w) all(w == 0))) == 0
#       is_not_dag = !is.dag(graph_from_adjacency_matrix(tree))
#       multi_parents = any(colSums(tree) > 1)
#
#       if(multi_root | no_root | is_not_dag | multi_parents) next;
#
#       # Pigeonhole constraint, it should have no violations of the principle
#       any_CCF_violation = sapply(rownames(tree), function(w)
#       {
#         ccf_node = df[w]
#
#         children = revolver:::children(tree, w)
#         ccf_children = sum(df[children])
#
#         ccf_node < ccf_children
#       })
#
#       if(any(any_CCF_violation)) next;
#
#       tree = revolver:::MatrixToDataFrame(tree)
#
#       span_trees = append(span_trees, list(tree))
#     }
#
#     all.trees = append(all.trees, list(span_trees))
#   }
#
#   names(all.trees) = colnames(m)
#
#   pio::pioStr("\nNumber of trees", paste(sapply(all.trees, length), collapse = ', '))
#
#   # print(lapply(all.trees, length))
#
#   all.trees
# }
#
# # data("Breast.fit")
# load("../test.revolver/TRACERx-release/TRACERx.cohort.RData", verbose = T)
#
# my.data = TRACERx.cohort$phylogenies$CRUK0001[[1]]$dataset
# sample.groups = TRACERx.cohort$phylogenies$CRUK0001[[1]]$samples_ID
# my.data2 = my.data
# my.data[, sample.groups] = my.data[, sample.groups] * 100
#
# clonal.cluster = '3'
#
# ct = revolver:::clusters.table(my.data2, sample.groups)
#
# my.obj = single_sample_CCF_pigeonhole_generator(ct[, sample.groups, drop = FALSE], 0.01)
#
# clonevol.obj = revolver:::useClonevo(my.data, sample.groups, clonal.cluster)
# my.data = revolver:::permuteClusterIds(my.data)
#
# capture.output({
#   clonevol.obj = revolver:::infer.clonal.models(
#     variants = my.data,
#     cluster.col.name = 'cluster',
#     # vaf.col.names = vaf.col.names,
#     ccf.col.names = sample.groups,
#     sample.names = sample.groups,
#     cancer.initiation.model = 'monoclonal',
#     # cancer.initiation.model = 'polyclonal',
#     # subclonal.test = 'bootstrap',
#     subclonal.test = 'none',
#     subclonal.test.model = 'non-parametric',
#     # subclonal.test.model = 'beta-binomial',
#     num.boots = 1000,
#     founding.cluster = clonal.cluster,
#     cluster.center = 'median',
#     ignore.clusters = NULL,
#     clone.colors = NULL,
#     min.cluster.vaf = 0.01,
#     # min probability that CCF(clone) is non-negative
#     sum.p = 0.05,
#     # alpha level in confidence interval estimate for CCF(clone)
#     alpha = 0.05,
#     verbose = F
#   )
# })
#
#
# # Queue functions
# Q = NULL
# pop = function(Q){ Q[-1] }
# push = function(Q,x){ append(Q, list(x)) }
# peek = function(Q){ Q[[1]] }
#
# # Data-query functions
# s = function(x){ df[df<x] }
#
# empty = matrix(0, nrow = length(df) + 1, ncol = length(df) + 1)
# colnames(empty) = rownames(empty) = c(names(df), "GL")
#
# empty["GL", rownames(ct)[ct$is.clonal]] = 1
#
# # Init is the clonal node
# Q = push(NULL, empty)
#
# # Loop untill Q is empty
# repeat{
#   if(length(Q) == 0) break;
#
#   # Get the head of the queue, it's a tree
#   head_node = peek(Q)
#
#   # Remove the head from the queue
#   Q = pop(Q)
#
#   # The head of the queue is a tree, we need to take all of
#   # its terminal nodes (leaves). We want only leaves with an
#   # incoming edge, which we get by looking at the transposed
#   leaves = setdiff(
#     revolver:::leaves(head_node),
#     revolver:::leaves(t(head_node)))
#
#   # Compute candidates nodes downstream of x
#   x = df[head_node]
#   candidates.x = s(x)
#
#   # Ensure that the candidates have not been processed
#   # visited = revolver:::MatrixToDataFrame(head_node)$to
#   # candidates.x =
#
#   # Filter those possible sets of children that
#   # violate the pigeonhole principle
#   children = lapply(1:length(candidates.x), function(w) {
#     com = combn(x = candidates.x, w, simplify = FALSE)
#     com = com[sapply(com, function(q) sum(q) <=  x)]
#
#     lapply(com, function(q) data.frame(from = names(x), to = names(q), stringsAsFactors = FALSE))
#     # com = lapply(com, function(q) data.frame(from = names(x), to = names(q), stringsAsFactors = FALSE))
#     # com[sapply(com, function(q) length(q) > 0)]
#   })
#
#   children = children[sapply(children, function(q) length(q) > 0)]
#   children = Reduce(append, children)
#
#   # Concatenate to the current head_node
#   head_node = revolver:::MatrixToDataFrame(head_node)
#   children = lapply(children, function(w) rbind(w, head_node))
#
#   # Convert to matrix
#   children = lapply(children, revolver:::DataFrameToMatrix)
#   for(ch in seq(children)) Q = push(Q, children[[ch]])
# }
#
# combn(candidates.x, 2, simplify = FALSE)
#
#
#
