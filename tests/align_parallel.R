library(parallel)
library(Biostrings)
library(DECIPHER)

# All cores

cores <- detectCores()

# Load raw alignment data

setwd("../data")

file_list <- list.files()
gisaid_fastas <- grep("gisaid_cov2020_sequences", file_list, value = TRUE)
unaligned <- readDNAStringSet(gisaid_fastas)

setwd("../tests")

# Define profile alignment function

align_profs <- function(x, y) {
  
  AlignProfiles(x, y, processors = NULL)
  
}

# Test 1 - Size 100 chunks 

tt1_1 <- Sys.time()

ua_subsets <- split(unaligned, ceiling(seq_along(unaligned)/100))

nogap_subsets <- mclapply(ua_subsets, function(x) RemoveGaps(x, removeGaps = "all", processors = NULL), mc.cores = cores)

align_subsets <- mclapply(nogap_subsets, function(x) AlignSeqs(x, iterations = 0, refinements = 0, processors = NULL), mc.cores = cores)

align_comp <- base::Reduce(align_profs, align_subsets)

tt1_2 <- Sys.time()

tt1_diff <- tt1_2 - tt1_1

print(tt1_diff)

print(aling_comp)

# Test 2 - Size 50 chunks 

tt2_1 <- Sys.time()

ua_subsets <- split(unaligned, ceiling(seq_along(unaligned)/50))

nogap_subsets <- mclapply(ua_subsets, function(x) RemoveGaps(x, removeGaps = "all", processors = NULL), mc.cores = cores)

align_subsets <- mclapply(nogap_subsets, function(x) AlignSeqs(x, iterations = 0, refinements = 0, processors = NULL), mc.cores = cores)

align_comp <- base::Reduce(align_profs, align_subsets)

tt2_2 <- Sys.time()

tt2_diff <- tt2_2 - tt2_1

print(tt2_diff)

print(aling_comp)

# Test 3 - Size 25 chunks 

tt3_1 <- Sys.time()

ua_subsets <- split(unaligned, ceiling(seq_along(unaligned)/25))

nogap_subsets <- mclapply(ua_subsets, function(x) RemoveGaps(x, removeGaps = "all", processors = NULL), mc.cores = cores)

align_subsets <- mclapply(nogap_subsets, function(x) AlignSeqs(x, iterations = 0, refinements = 0, processors = NULL), mc.cores = cores)

align_comp <- base::Reduce(align_profs, align_subsets)

tt3_2 <- Sys.time()

tt3_diff <- tt3_2 - tt3_1

print(tt3_diff)

print(aling_comp)

# Test 4 - Size 10 chunks 

tt4_1 <- Sys.time()

ua_subsets <- split(unaligned, ceiling(seq_along(unaligned)/10))

nogap_subsets <- mclapply(ua_subsets, function(x) RemoveGaps(x, removeGaps = "all", processors = NULL), mc.cores = cores)

align_subsets <- mclapply(nogap_subsets, function(x) AlignSeqs(x, iterations = 0, refinements = 0, processors = NULL), mc.cores = cores)

align_comp <- base::Reduce(align_profs, align_subsets)

tt4_2 <- Sys.time()

tt4_diff <- tt4_2 - tt4_1

print(tt4_diff)

print(aling_comp)
