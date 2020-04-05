library(data.table)
library(plyr)
library(RColorBrewer)

# Load testing data 

loadRData <- function(fileName){
load(fileName)
get(ls()[ls() != "fileName"])
}

setwd("../data")

file_list <- list.files() 

pre_aligned_filtered <- loadRData(grep("dec_aligned_filtered", file_list, value = TRUE))
pre_meta <- loadRData(grep("covid_filtered_meta", file_list, value = TRUE))
pre_dist <- loadRData(grep("dec_fasta_dist", file_list, value = TRUE))
meta <- as.data.frame(pre_meta[,c("Accession", "Region", "Geo_Location", "Datetime")])
vars_freq <- loadRData(grep("var_freq_sub*", file_list, value = TRUE))

setwd("..")

# Load color palettes    

kev_palette <- c(
"dodgerblue2", "#E31A1C",
"green4",
"#6A3D9A", 
"#FF7F00", 
"black", "gold1",
"skyblue2", "#FB9A99", 
"palegreen2",
"#CAB2D6", 
"#FDBF6F", 
"gray70", "khaki2",
"maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
"darkturquoise", "green1", "yellow4", "yellow3",
"darkorange4", "brown"
)

qual_palettes = brewer.pal.info[brewer.pal.info$category == "qual", ]
qual_vector = unlist(mapply(brewer.pal, qual_palettes$maxcolors, rownames(qual_palettes)))

# Source global

source("global.R")

setwd("data")

# umap

umap_preloaded <- umap_process(pre_dist, meta)

save(umap_preloaded, file = "umap_preloaded.RData")

# mst

mst_preloaded <- mst_graph(pre_dist, meta)

save(mst_preloaded, file = "mst_preloaded.RData")

# snps
  
snps_preloaded <- snps_get(pre_aligned_filtered, meta, vars_freq)

save(snps_preloaded, file = "snps_preloaded.RData")


