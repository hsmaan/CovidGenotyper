library(uwot)

umap_process <- function(covid_dist) {
  
  covid_dist <- dist(covid_dist)
  set.seed(2020)
  covid_umap <- uwot::umap(covid_dist, init = "spectral", metric = "cosine", n_neighbors = 50, min_dist = 0.001, spread = 30, local_connectivity = 10)
  covid_umap_df <- as.data.frame(covid_umap)
  return(covid_umap_df)
  
}
