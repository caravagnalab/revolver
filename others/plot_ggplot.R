install.packages('tidygraph', dependencies = T)
library(tidygraph)
library(ggraph)


data('Breast.fit', package = 'revolver')

tree = Breast.fit$phylogenies[[20]][[1]]
nodes.annot = tree$CCF[, 1:3]
nodes.annot$name = rownames(nodes.annot)

revolver:::plot.rev_phylo(tree)

tidy = as_tbl_graph(tree$adj_mat)

tidy = tidy %>% activate(nodes) %>%
  left_join(nodes.annot, by = 'name') %>%
  mutate(
    is.driver = ifelse(name == 'GL', FALSE, is.driver),
    is.clonal = ifelse(name == 'GL', FALSE, is.clonal),
    nMuts = ifelse(name == 'GL', 0, nMuts),
  )

colors = rep('gainsboro', )

  mobster:::scols(rownames(tree$CCF)[tree$CCF$is.driver])


edges_colors = revolver:::MatrixToDataFrame(tree$adj_mat)
transfers = tree$transfer$clones

edColor.path = NULL

for (i in 1:nrow(transfers))
{
  # Get the path
  p = revolver:::find.path(tree$adj_mat, from = transfers[i, 1], to = transfers[i, 2])
  
}
  

transfers = as_tibble(transfers) %>% filter(from != 'GL')

transfers$from = as.integer(transfers$from)
transfers$to = as.integer(transfers$to)


  scale_edge_colour_manual(values = c(ally = "#22B022",
                                      opponent = "#A4AAF6")) 
  # tidy %>% activate(nodes) %>%
  # mutate(Species = ifelse(leaf, as.character(iris$Species)[label], NA)) %>% 
  #   activate(edges) %>% 
  #   mutate(to_setose = .N()$Species[to] == 'setosa')
  
tidy %>% activate(edges) %>%
  left_join(tree$transfer$clones, by = c('from', 'to')) %>%
  mutate(
    is.driver = ifelse(name == 'GL', FALSE, is.driver),
    is.clonal = ifelse(name == 'GL', FALSE, is.clonal),
    nMuts = ifelse(name == 'GL', 0, nMuts),
  )
  
tidy %>% 
  # mutate(leaf = node_is_leaf(), root = node_is_root()) %>% 
  ggraph(layout = 'tree') +
  # geom_edge_diagonal() +
  # geom_node_point(aes(filter = is.driver), colour = 'forestgreen', size = 10) +
  geom_node_point(aes(filter = !is.driver, size = nMuts), colour = 'gainsboro') +
  # geom_node_point(size = 6, colour = 'gainsboro') +
  geom_node_point(aes(filter = is.driver, colour = name, size = nMuts)) +
  geom_node_point(aes(filter = (name == 'GL')), colour = 'white', size = 10) +
  geom_node_text(aes(label = name), colour = 'black', vjust = 0.4) + 
  geom_node_text(aes(filter = (name == 'GL'), label = name), colour = 'steelblue', vjust = 0.4) + 
  scale_size_continuous(range = c(6, 10)) +
  geom_edge_link(
    # aes(colour = type),
    arrow = arrow(length = unit(1.5, "mm")),
    start_cap = circle(3, "mm"),
    end_cap = circle(3, "mm")
  ) + 
  guides(color = FALSE) + 
  scale_color_brewer(palette = 'Set1') + 
  labs(size = "nMuts") +
  # geom_node_point(aes(filter = !is.driver), colour = 'gray', size = 5) +
  # geom_node_point(aes(filter = root), colour = 'firebrick', size = 10) +
  theme_graph() +
  theme(legend.position="bottom")



tidy %>% 
  mutate(leaf = node_is_leaf(), root = node_is_root()) %>% 
  ggraph(layout = 'tree') +
  geom_edge_diagonal2(n = 2)+
  geom_node_point(aes(filter = leaf), colour = 'forestgreen', size = 10) +
  geom_node_point(aes(filter = !leaf), colour = 'gray', size = 5) +
  geom_node_point(aes(filter = root), colour = 'firebrick', size = 10) +
  theme_graph()
