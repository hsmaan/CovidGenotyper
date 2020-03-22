library(stringr)

multi_dim_scale <- function(covid_dist) {
  
  cmd <- cmdscale(d = covid_dist, k = 2)
  acc_names = rownames(cmd)
  cmd <- data.frame("Accession" = acc_names, "MDS_1" = cmd[,1], "MDS_2" = cmd[,2])
  return(cmd)
  
}
