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

# Get complete results from dist.dna call

attach(loadNamespace("ape"), name = "ape_ns")

dist_dna_mod <- function (x, model = "K80", variance = FALSE, gamma = FALSE, 
                          pairwise.deletion = FALSE, base.freq = NULL, as.matrix = FALSE) 
{
  MODELS <- c("RAW", "JC69", "K80", "F81", "K81", "F84", "T92", 
              "TN93", "GG95", "LOGDET", "BH87", "PARALIN", "N", "TS", 
              "TV", "INDEL", "INDELBLOCK")
  imod <- pmatch(toupper(model), MODELS)
  if (is.na(imod)) 
    stop(paste("'model' must be one of:", paste("\"", MODELS, 
                                                "\"", sep = "", collapse = " ")))
  if (imod == 11 && variance) {
    warning("computing variance not available for model BH87")
    variance <- FALSE
  }
  if (gamma && imod %in% c(1, 5:7, 9:17)) {
    warning(paste("gamma-correction not available for model", 
                  model))
    gamma <- FALSE
  }
  if (is.list(x)) 
    x <- as.matrix(x)
  nms <- dimnames(x)[[1]]
  n <- dim(x)
  s <- n[2]
  n <- n[1]
  if (s * n > 2^31 - 1) 
    stop("dist.dna() cannot handle more than 2^31 - 1 bases")
  if (imod %in% c(4, 6:8)) {
    BF <- if (is.null(base.freq)) 
      base.freq(x)
    else base.freq
  }
  else BF <- 0
  if (imod %in% 16:17) 
    pairwise.deletion <- TRUE
  if (!pairwise.deletion) {
    keep <- .C(GlobalDeletionDNA, x, n, s, rep(1L, s))[[4]]
    x <- x[, as.logical(keep)]
    s <- dim(x)[2]
  }
  Ndist <- if (imod == 11) 
    n * n
  else n * (n - 1)/2
  var <- if (variance) 
    double(Ndist)
  else 0
  if (!gamma) 
    gamma <- alpha <- 0
  else {
    alpha <- gamma
    gamma <- 1
  }
  d <- .C(dist_dna, x, as.integer(n), as.integer(s), imod, 
          double(Ndist), BF, as.integer(pairwise.deletion), as.integer(variance), 
          var, as.integer(gamma), as.double(alpha), NAOK = TRUE)
  return(d)
  # if (variance) 
  #   var <- d[[9]]
  # d <- d[[5]]
  # if (imod == 11) {
  #   dim(d) <- c(n, n)
  #   dimnames(d) <- list(nms, nms)
  # }
  # else {
  #   attr(d, "Size") <- n
  #   attr(d, "Labels") <- nms
  #   attr(d, "Diag") <- attr(d, "Upper") <- FALSE
  #   attr(d, "call") <- match.call()
  #   attr(d, "method") <- model
  #   class(d) <- "dist"
  #   if (as.matrix) 
  #     d <- as.matrix(d)
  # }
  # if (variance) 
  #   attr(d, "variance") <- var
  # d
}

# Test for subset of pre_aligned

pre_aligned_sub <- as.DNAbin(pre_aligned_filtered[1:50])

test1_dist <- dist.dna(pre_aligned_sub, model = "K80", as.matrix = TRUE, pairwise.deletion = FALSE)

test2_dist <- dist_dna_mod(pre_aligned_sub, model = "K80", as.matrix = TRUE, pairwise.deletion = FALSE)

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
  d2 <- DNAStringSet(d2)
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

align_mat_num <- alignment_mat

align_mat_num[align_mat_num == "A"] <- 1
align_mat_num[align_mat_num == "G"] <- 2
align_mat_num[align_mat_num == "T"] <- 3
align_mat_num[align_mat_num == "C"] <- 4
align_mat_num[align_mat_num == "N"] <- 5
align_mat_num[align_mat_num == "-"] <- 0

align_mat_sub <- align_mat_num[1:200,]
align_mat_sub_num <- apply(align_mat_sub, 1, as.numeric)

# Test paralleldist implementation

profvis({
  test <- dist.dna(as.DNAbin(pre_aligned_filtered), model = "raw", as.matrix = TRUE)
})
# 
min(test_combine)

# Test distance update 

align <- pre_aligned_filtered
# new_length <- length(readDNAStringSet(fasta))
new_length <- 4
align_new <- align[(length(align) - new_length + 1):length(align)]
align_minus <- align[1:(length(align) - new_length)]

align_new_list <- as.list(align_new)


split_test <- lapply(align_new, function(x) lapply(align_minus, function(y) split_dist(x, y)))

split_test_reduced <- lapply(split_test, function(x) Reduce(cbind, x))

split_test_reduced_mins <- lapply(split_test_reduced, min)
split_test_reduced_min_id <- lapply(split_test_reduced, which.min)
split_test_reduced_mins <- lapply(split_test_reduced, min)
split_test_reduced_names <- as.list(names(split_test_reduced_mins))

align_minus_up <- align_minus
pre_dist_minus <- pre_dist[1:(nrow(pre_dist) - new_length), 1:(ncol(pre_dist) - new_length)]

align_update <- function(min_id, name) {
  
  dist_min_row <- pre_dist_minus[min_id[1],]
  return(list(dist_min_row, name))
  
}

align_bind <- function(full_dist, row_list) {
  
  new_row <- row_list[[1]]
  name <- row_list[[2]]
  new_col <- c(new_row, 0)
  dist_new <- rbind(full_dist, new_row)
  rownames(dist_new)[nrow(dist_new)] <- name
  dist_new <- cbind(dist_new, new_col)
  colnames(dist_new)[ncol(dist_new)] <- name
  return(dist_new)
}

test <- mapply(align_update, min_id = split_test_reduced_min_id, name = split_test_reduced_names, SIMPLIFY = FALSE)

test2 <- c(list(pre_dist_minus), test)

test3 <- Reduce(align_bind, test2)



align_update(split_test_reduced_min_id[[1]], split_test_reduced_names[[1]])

if (max(unlist(split_test_reduced_mins)) > 1e-4) {
  print("test")
} else {
  
  
  
}



split_dist <- function(d1, d2) {
  
  dna_bound <- c(d1, d2)
  bound_dist <- dist.dna(dna_bound, as.matrix = TRUE)
  bound_df <- as.data.frame(bound_dist[2])
  colnames(bound_df) <- names(dna_bound)[1]
  rownames(bound_df) <- names(dna_bound)[2]
  return(bound_df)
  
}
