library(stringr)

multi_dim_scale <- function(covid_dist) {
  
  cmd <- cmdscale(d = covid_dist, k = 2)
  acc_names = (str_split_fixed(rownames(cmd), fixed("."), 2)[,1])
  cmd <- data.frame("Accession" = acc_names, "MDS_1" = cmd[,1], "MDS_2" = cmd[,2])
  return(cmd)
  
}
