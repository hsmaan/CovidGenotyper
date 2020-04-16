library(data.table)
library(plyr)
library(RColorBrewer)
library(profvis)
library(parallel)

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

# Source global

source("R/global.R")

setwd("data")

# Load ape 0.1 function

dist.dna.Kimura <- function(x, y, variance = FALSE, gamma = NULL) {
  L <- length(x)
  d <- x != y
  Nd <- sum(d)
  pw.diff <- cbind(x[d], y[d])
  PuPy <- ifelse(pw.diff == "a" | pw.diff == "g", "R", "Y")
  Nv <- sum(PuPy[, 1] != PuPy[, 2])
  Ns <- Nd - Nv
  P <- Ns / L
  Q <- Nv / L
  a1 <- 1 - 2 * P - Q
  a2 <- 1 - 2 * Q
  if (is.null(gamma)) D <- -0.5 * log(a1 * sqrt(a2))
  else {
    b <- -1 / gamma
    D <- gamma * (a1^b + 0.5 * a2^b - 1.5) / 2
  }
  if (variance) {
    if (is.null(gamma)) {
      c1 <- 1 / a1
      c2 <- 1 / a2
      c3 <- (c1 + c2) / 2
    }
    else {
      b <- -(1 / gamma + 1)
      c1 <- a1^b
      c2 <- a2^b
      c3 <- (c1 + c2) / 2            
    }
    var.D <- (c1^2 * P + c3^2 * Q - (c1 * P + c3 * Q)^2) / L
    return(c(D, var.D))
  }
  else return(D)
}

# Subset alignments

alignment_minus <- pre_aligned_filtered[-length(pre_aligned_filtered)]
alignment_last <- pre_aligned_filtered[length(pre_aligned_filtered)]

# Convert to matrix

align_minus_mat <- as.matrix(alignment_minus)
align_last_vec <- as.vector(alignment_last)

# Get distance 

split_dist <- function(d1, d2) {
  
  d1 <- DNAStringSet(d1)
  dna_bound <- c(d1, d2)
  bound_dist <- dist_get(dna_bound)
  bound_df <- as.data.frame(bound_dist[2])
  colnames(bound_df) <- names(dna_bound)[1]
  rownames(bound_df) <- names(dna_bound)[2]
  return(bound_df)
  
}

profvis({
  split_test <- lapply(alignment_minus, function(x) split_dist(x, alignment_last))
  test_combine <- as.vector(unlist(Reduce(cbind, split_test)))
})

# Get numeric vector from matrix

alignment_mat <- as.matrix(pre_aligned_filtered)

alignment_mat_num <- sapply(alignment_mat, function(x) plyr::mapvalues(x, from = c("A", "T", "C", "G", "N"), to = c(1, 2, 3, 4, 5)))

