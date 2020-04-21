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
pre_umap <- loadRData(grep("umap_preloaded", file_list, value = TRUE))
meta <- as.data.frame(pre_meta[,c("Accession", "Region", "Geo_Location", "Datetime")])
vars_freq <- loadRData(grep("var_freq_sub*", file_list, value = TRUE))

setwd("..")

# Source global

source("R/global.R")

setwd("data")

# Diff between meta and new

acc_names = rownames(covid_dist)

umap_update <- function(min_id, umap_pre, number) {
  
  umap_min_row <- umap_pre[min_id[1],]
  umap_new <- rbind(umap_pre, umap_min_row)
  umap_new$Accession[nrow(umap_new)] <- new_acc_names[number]
  umap_ret <<- umap_new
  
}

pre_aligned_filtered <- pre_aligned_filtered[order(names(pre_aligned_filtered))]
umap_0 <- pre_umap
umap_0 <- umap_0[order(umap_0$Accession),]
# new_length <- length(readDNAStringSet(fasta))
new_length <- 2
acc_names <- rownames(pre_dist)
new_acc_names <- acc_names[(length(acc_names) - new_length):length(acc_names)] 
new_length_num <- as.list(seq(new_length))
align_new <- pre_aligned_filtered[(length(pre_aligned_filtered) - new_length + 1):length(pre_aligned_filtered)]
align_minus <- pre_aligned_filtered[1:(length(pre_aligned_filtered) - new_length)]
align_new_list <- as.list(align_new)

split_test <- mclapply(align_new, function(x) mclapply(align_minus, function(y) split_dist(x, y), mc.cores = cores), mc.cores = cores)
split_test_reduced <- mclapply(split_test, function(x) Reduce(cbind, x), mc.cores = cores)
split_test_reduced_mins <- mclapply(split_test_reduced, min, mc.cores = cores)
split_test_reduced_min_id <- mclapply(split_test_reduced, which.min, mc.cores = cores)
split_test_reduced_names <- as.list(names(split_test_reduced_mins))

umap_ret <- umap_minus
umap_update <- function(min_id, umap_pre, number) {
  umap_min_row <- umap_pre[min_id[1],]
  umap_new <- rbind(umap_pre, umap_min_row)
  umap_new$Accession[nrow(umap_new)] <- new_acc_names[number]
  umap_ret <<- umap_new
}
mcmapply(umap_update, min_id = split_test_reduced_min_id, number = new_length_num, MoreArgs = list(umap_pre = umap_ret), mc.cores = cores)
return(umap_ret)

umap_update(split_test_reduced_min_id[[1]], umap_pre = umap_minus, number = new_length_num[[1]])

umap_min_row <- umap_minus[split_test_reduced_min_id[[1]][1],]
umap_new <- rbind(umap_minus, umap_min_row)
umap_new$Accession[nrow(umap_new)] <- new_acc_names[1]
umap_ret <<- umap_new

umap_process_heur <- function(align, fasta, new_dist, new_meta, old_umap) {
  
  umap_0 <- old_umap
  umap_0$Accession <- as.character(umap_0$Accession)
  new_length <- length(readDNAStringSet(fasta))
  acc_names <- new_meta$Accession
  new_acc_names <- acc_names[(length(acc_names) - new_length + 1):length(acc_names)]
  meta_adds <- new_meta[(length(acc_names) - new_length + 1):length(acc_names),]
  new_length_num <- as.list(seq(new_length))
  align_new <- align[(length(align) - new_length + 1):length(align)]
  align_minus <- align[1:(length(align) - new_length)]
  align_new_list <- as.list(align_new)
  
  split_test <- mclapply(align_new, function(x) mclapply(align_minus, function(y) split_dist(x, y), mc.cores = cores), mc.cores = cores)
  split_test_reduced <- mclapply(split_test, function(x) Reduce(cbind, x), mc.cores = cores)
  split_test_reduced_mins <- mclapply(split_test_reduced, min, mc.cores = cores)
  split_test_reduced_min_id <- mclapply(split_test_reduced, which.min, mc.cores = cores)
  split_test_reduced_names <- as.list(names(split_test_reduced_mins))
  
  if (max(unlist(split_test_reduced_mins)) > 1e-4) {
    umap_ret <- umap_process(new_dist, meta_data)
    return(umap_ret)
  } else {
    umap_ret <- umap_0
    umap_update <- function(min_id, number) {
      umap_min_row <- umap_ret[min_id[1],]
      umap_min_row <- umap_min_row[,1:3]
      umap_min_row$Accession[1] <- new_acc_names[number]
      min_row_meta <- meta_adds[number,]
      umap_min_final <- merge(umap_min_row, min_row_meta)
      colnames(umap_min_final) <- colnames(umap_ret)
      return(umap_min_final)
    }
    umap_updates <- mcmapply(umap_update, min_id = split_test_reduced_min_id, number = new_length_num, mc.cores = cores, SIMPLIFY = FALSE)
    umap_updates <- c(list(umap_ret), umap_updates)
    umap_ret <- Reduce(rbind, umap_updates)
    return(umap_ret)
  }
}

# Test (old align (-5), old_umap (-5), new meta and dist)

align_minus <- pre_aligned_filtered[1:(length(pre_aligned_filtered) - 2)]

umap_minus <- pre_umap[1:(nrow(pre_umap) - 2),]

test <- umap_process_heur(pre_aligned_filtered, "test", pre_dist, meta, umap_minus)

      dist_get_heur <- function(align, fasta, dist) {
  
  dist_ret0 <- dist
  new_length <- fasta
  align_new <- align[(length(align) - new_length + 1):length(align)]
  align_minus <- align[1:(length(align) - new_length)]
  align_new_list <- as.list(align_new)
  
  split_test <- mclapply(align_new, function(x) mclapply(align_minus, function(y) split_dist(x, y), mc.cores = cores), mc.cores = cores)
  split_test_reduced <- mclapply(split_test, function(x) Reduce(cbind, x), mc.cores = cores)
  split_test_reduced_mins <- mclapply(split_test_reduced, min, mc.cores = cores)
  split_test_reduced_min_id <- mclapply(split_test_reduced, which.min, mc.cores = cores)
  split_test_reduced_names <- as.list(names(split_test_reduced_mins))
  
  if (max(unlist(split_test_reduced_mins)) > 1e-4) {
    dist_ret <- dist_get(align)
    return(dist_ret)
  } else {
    dist_ret <- dist_ret0
    align_update <- function(min_id, name, dist) {
      dist_min_row <- dist[min_id[1],]
      dist_min_col <- c(dist_min_row, 0)
      dist_new <- rbind(dist, dist_min_row)
      rownames(dist_new)[nrow(dist_new)] <- name
      dist_new <- cbind(dist_new, dist_min_col)
      colnames(dist_new)[ncol(dist_new)] <- name
      dist_ret <<- dist_new
    }
    mcmapply(align_update, min_id = split_test_reduced_min_id, name = split_test_reduced_names, MoreArgs = list(dist = dist_ret), mc.cores = cores)
    return(dist_ret)
  }
}

test <- dist_get_heur(pre_aligned_filtered, 2, pre_dist_minus)
