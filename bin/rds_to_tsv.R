library(data.table)

setwd("../data")

# Load files

files <- list.files()

loadRData <- function(fileName){
	  load(fileName)
  get(ls()[ls() != "fileName"])
}

dist_rds <- loadRData(grep("dec_fasta_dist*", files, value = TRUE))
meta_rds <- loadRData(grep("covid_filtered_meta_*", files, value = TRUE))
umap_rds <- loadRData(grep("umap_preloaded*", files, value = TRUE))


# Reformat into tsvs

dist_dt <- as.data.table(dist_rds, keep.rownames = TRUE)

meta_dt <- as.data.table(meta_rds)

umap_dt <- as.data.table(umap_rds)

# Remove old files

file.remove(grep("signal_dist_", files, value = TRUE))
file.remove(grep("signal_meta_", files, value = TRUE))
file.remove(grep("signal_umap_", files, value = TRUE))

# Output new files

fwrite(dist_dt, file = paste("signal_dist_", Sys.Date(), ".tsv", sep = ""), quote = FALSE, col.names = TRUE, row.names = FALSE)
fwrite(meta_dt, file = paste("signal_meta_", Sys.Date(), ".tsv", sep = ""), quote = FALSE, col.names = TRUE, row.names = FALSE)
fwrite(umap_dt, file = paste("signal_umap_", Sys.Date(), ".tsv", sep = ""), quote = FALSE, col.names = TRUE, row.names = FALSE)


