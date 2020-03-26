library(igraph)

mst_graph <- function(covid_dist, meta_data, vertex_cols) {
  
  g <- graph.adjacency(covid_dist, mode = "undirected", weighted = TRUE, diag = FALSE)
  g_mst <- mst(g)
  acc_ordering <- match(meta_data[,1], names(V(g_mst)[[]]))
  meta_ordered <- meta_data[order(acc_ordering), ]
  meta_colors <- (vertex_cols[1:length(unique(meta_ordered[,2]))])[factor(meta_ordered[,2])]
  V(g_mst)$color <- meta_colors
  return(g_mst)
  
}
