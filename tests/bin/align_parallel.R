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

ua_subsets <- split(unaligned, ceiling(seq_along(unaligned)/50))

nogap_subsets <- mclapply(ua_subsets, function(x) RemoveGaps(x, removeGaps = "all", processors = 2), mc.cores = cores)

align_subsets <- mclapply(nogap_subsets, function(x) AlignSeqs(x, iterations = 0, refinements = 0, processors = 2), mc.cores = cores)

align_comp <- base::Reduce(align_profs, align_subsets)

tt1_2 <- Sys.time()

tt1_diff <- tt1_2 - tt1_1

print(tt1_diff)

print(align_comp)

# Test 2 - Size 50 chunks 

#tt2_1 <- Sys.time()

#ua_subsets_2 <- split(unaligned, ceiling(seq_along(unaligned)/50))

#nogap_subsets_2 <- mclapply(ua_subsets_2, function(x) RemoveGaps(x, removeGaps = "all", processors = 2), mc.cores = cores)

#align_subsets_2 <- mclapply(nogap_subsets_2, function(x) AlignSeqs(x, iterations = 0, refinements = 0, processors = 2), mc.cores = cores)

#align_comp_2 <- base::Reduce(align_profs, align_subsets_2)

#tt2_2 <- Sys.time()

#tt2_diff <- tt2_2 - tt2_1

#print(tt2_diff)

#print(align_comp_2)

# Test 3 - Size 25 chunks 

#tt3_1 <- Sys.time()

#ua_subsets_3 <- split(unaligned, ceiling(seq_along(unaligned)/25))

#nogap_subsets_3 <- mclapply(ua_subsets_3, function(x) RemoveGaps(x, removeGaps = "all", processors = 2), mc.cores = cores)

#align_subsets_3 <- mclapply(nogap_subsets_3, function(x) AlignSeqs(x, iterations = 0, refinements = 0, processors = 2), mc.cores = cores)

#align_comp_3 <- base::Reduce(align_profs, align_subsets_3)

#tt3_2 <- Sys.time()

#tt3_diff <- tt3_2 - tt3_1

#print(tt3_diff)

#print(align_comp_3)

# Test 4 - Size 10 chunks 

#tt4_1 <- Sys.time()

#ua_subsets_4 <- split(unaligned, ceiling(seq_along(unaligned)/10))

#nogap_subsets_4 <- mclapply(ua_subsets_4, function(x) RemoveGaps(x, removeGaps = "all", processors = 2), mc.cores = cores)

#align_subsets_4 <- mclapply(nogap_subsets_4, function(x) AlignSeqs(x, iterations = 0, refinements = 0, processors = 2), mc.cores = cores)

#align_comp_4 <- base::Reduce(align_profs, align_subsets_4)

#tt4_2 <- Sys.time()

#tt4_diff <- tt4_2 - tt4_1

#print(tt4_diff)

#print(align_comp_4)
