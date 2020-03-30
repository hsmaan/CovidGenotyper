library(data.table)
library(stringr)

setwd("../data")

all_files <- list.files()

meta_file <- grep("gisaid_metadata", all_files, value = TRUE)

metadata <- fread(meta_file, stringsAsFactors = FALSE)

meta_df <- as.data.frame(metadata)

meta_sub <- meta_df[,c("date", "gisaid_epi_isl", "region", "country", "length", "age", "sex")]

colnames(meta_sub) <- c("Date", "Accession", "Region", "Geo_Location", "Genome_Length", "Age", "Sex")

save(meta_sub, file = paste("covid_meta_", Sys.Date(), ".RData", sep = ""))
