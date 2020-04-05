library(data.table)
library(stringr)

setwd("../data")

all_files <- list.files()

meta_file <- grep("gisaid_metadata", all_files, value = TRUE)
  
metadata <- fread(meta_file, stringsAsFactors = FALSE)

meta_df <- as.data.frame(metadata)

meta_sub <- meta_df[,c("date", "gisaid_epi_isl", "region", "country", "length", "age", "sex")]

colnames(meta_sub) <- c("Date", "Accession", "Region", "Geo_Location", "Genome_Length", "Age", "Sex")

'%ni%' <- Negate('%in%')

meta_sub <- meta_sub[meta_sub$Region != " ", ]

meta_sub <- meta_sub[meta_sub$Geo_Location != " ", ]

meta_sub <- meta_sub[(str_length(meta_sub$Date) > 9),]

meta_sub <- meta_sub[(rownames(meta_sub) %ni% (grep("XX", meta_sub$Date))),]

sample_time <- meta_sub$Date
sample_time <- as.Date(sample_time, format = "%Y-%m-%d")
sample_time <- unclass(sample_time)
orig_time <- as.Date("2019-12-01", format = "%Y-%m-%d")
orig_time <- unclass(orig_time)
sample_time <- sample_time - orig_time

current_orig <- (unclass(Sys.Date())) - (unclass(as.Date("2019-12-01", format = "%Y-%m-%d")))

meta_sub$Datetime <- sample_time
meta_sub <- meta_sub[(meta_sub$Datetime >= 0), ]
meta_sub <- meta_sub[(meta_sub$Datetime <= current_orig), ]

save(meta_sub, file = paste("covid_meta_", Sys.Date(), ".RData", sep = ""))
