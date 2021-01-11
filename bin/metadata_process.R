library(data.table)
library(stringr)

setwd("../data")

# Load metadata

all_files <- list.files()
meta_file <- grep("gisaid_metadata", all_files, value = TRUE)
metadata <- fread(meta_file, stringsAsFactors = FALSE)

# Subset appropriate columns

meta_df <- as.data.frame(metadata)
meta_sub <- meta_df[,c(
    "date",
    "gisaid_epi_isl",
    "strain",
    "region", 
    "country", 
    "length", 
    "age", 
    "sex", 
    "country_exposure", 
    "pangolin_lineage"
)]
colnames(meta_sub) <- c(
    "Date", 
    "Accession",
    "Name", 
    "Region", 
    "Geo_Location", 
    "Genome_Length", 
    "Age", 
    "Sex", 
    "Country_Exposure",
    "Pangolin_Lineage"
)

# Define negate 

'%ni%' <- Negate('%in%')

# Format columns 

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

meta_sub$Country_Exposure <- ifelse((meta_sub$Geo_Location == meta_sub$Country_Exposure), "Not available", meta_sub$Country_Exposure)

# Subset for dates after Sept 1, 2020
cutoff_date <- (unclass(as.Date("2020-09-01", format = "%Y-%m-%d")) - (unclass(as.Date("2019-12-01", format = "%Y-%m-%d"))))
meta_sub <- meta_sub[which(meta_sub$Datetime > cutoff_date), ]

# Subset ack columns

meta_ack <- meta_df[,c("gisaid_epi_isl", "originating_lab", "submitting_lab", "authors", "date_submitted")]
colnames(meta_ack) <- c("Accession", "Originating lab", "Submitting lab", "Authors", "Date submitted")

# Save files

file.remove(grep("covid_meta_*", all_files, value = TRUE))
file.remove(grep("../ack/gisaid_acknowledgements_*", all_files, value = TRUE))
save(meta_sub, file = paste("covid_meta_", Sys.Date(), ".RData", sep = ""))
fwrite(meta_ack, file = paste("../ack/gisaid_acknowledgements_", Sys.Date(), ".csv", sep = ""), sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)
