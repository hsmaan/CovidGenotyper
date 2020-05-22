library(Biostrings)
library(stringr)
library(DECIPHER)
library(ape)
library(data.table)
library(uwot)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)

# Functions from CovidGenotyper/R/global.R  

align_get <- function(fasta, align) {
  
  covid_seq <- fasta
  if (length(covid_seq) > 10) {
    stop("Too many sequences in fasta file, only up to 10 viral genomes allowed at a time.")
  }
  if (length(covid_seq) < 1) {
    stop("No sequences found in fasta file.")
  }
  if (length(which((mclapply(covid_seq, length, mc.cores = cores)) < 29000)) > 0) {
    stop("One or more sequences not complete (length < 29000 nucleotides).")
  }
  covid_align <- align
  covid_seq <- RemoveGaps(covid_seq, removeGaps = "all", processors = NULL)
  fasta_final <- AlignProfiles(covid_align, covid_seq, processors = NULL)
  return(fasta_final)
  
}

dist_get <- function(align) {
  
  mask_sites <- c(187, 1059, 2094, 3037, 3130, 6990, 8022, 10323, 10741, 11074, 13408, 14786, 19684, 20148, 21137, 24034, 24378, 25563, 26144, 26461, 26681, 28077, 28826, 28854, 29700, 4050, 13402, 11083, 15324, 21575)
  align_mat <- as.matrix(align)
  align_mat_sub <- align_mat[, -mask_sites]
  align_mat_bin <- as.DNAbin(align_mat_sub)
  align_masked <- align_mat_bin %>% as.list %>% as.character %>% lapply(., paste0, collapse = "") %>% unlist %>% DNAStringSet
  align_trim <- subseq(align_masked, start = 265, end = 29674)
  dec_dist <- dist.dna(as.DNAbin(align_trim), model = "K80", as.matrix = TRUE, pairwise.deletion = FALSE)
  colnames(dec_dist) <- (str_split_fixed(colnames(dec_dist), fixed("."), 2)[,1])
  rownames(dec_dist) <- (str_split_fixed(rownames(dec_dist), fixed("."), 2)[,1])
  return(dec_dist)
  
}

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

align_update <- function(dist_mat, min_id, name) {
  
  dist_min_row <- dist_mat[min_id[1],]
  return(list(dist_min_row, name))
  
}

dist_get_heur <- function(align, fasta, dist) {
  
  dist_ret0 <- dist
  new_length <- length(fasta)
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
    update_rows <- mcmapply(align_update, min_id = split_test_reduced_min_id, name = split_test_reduced_names, MoreArgs = list(dist_mat = dist_ret), SIMPLIFY = FALSE, mc.cores = cores)
    dist_list <- c(list(dist_ret), update_rows)
    dist_concat <- Reduce(align_bind, dist_list)
    return(dist_concat)
  }
}

umap_process <- function(covid_dist, meta_df) {
  
  acc_names = rownames(covid_dist)
  covid_dist <- dist(covid_dist)
  set.seed(2020)
  covid_umap <- uwot::umap(covid_dist, init = "spectral", metric = "cosine", n_neighbors = 50, min_dist = 0.001, spread = 40, local_connectivity = 10, n_threads = cores*2)
  covid_umap_df <- as.data.frame(covid_umap)
  umap_df_final <- data.frame("Accession" = acc_names, "UMAP_1" = covid_umap_df[,1], "UMAP_2" = covid_umap_df[,2])
  umap_df_final <- merge(umap_df_final, meta_df)
  colnames(umap_df_final) <- c("Accession", "UMAP_1", "UMAP_2", "Region", "Country", "Date")
  return(umap_df_final)
  
}

umap_process_heur <- function(align, fasta, new_dist, new_meta, old_umap) {
  
  umap_0 <- old_umap
  umap_0$Accession <- as.character(umap_0$Accession)
  new_length <- length(fasta)
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

# Load color palettes    

kev_palette <- c("dodgerblue2", 
		 "#E31A1C",
		 "green4",
		 "#6A3D9A", 
	         "#FF7F00", 
		 "black",
		 "gold1",
		 "skyblue2",
		 "#FB9A99", 
		 "palegreen2",
		 "#CAB2D6", 
		 "#FDBF6F", 
		 "gray70", 
		 "khaki2",
		 "maroon",
		 "orchid1",
		 "deeppink1",
		 "blue1",
		 "steelblue4",
	         "darkturquoise",
	       	 "green1",
		 "yellow4",
		 "yellow3",
		 "darkorange4",
		 "brown")

qual_palettes = brewer.pal.info[brewer.pal.info$category == "qual", ]
qual_vector = unlist(mapply(brewer.pal, qual_palettes$maxcolors, rownames(qual_palettes)))

# Read args

args <- commandArgs(trailingOnly = TRUE)

fasta_file <- as.character(args[1])
align_file <- as.character(args[2])
dist_file <- as.character(args[3])
meta_file <- as.character(args[4])
umap_file <- as.character(args[5])
sample_hist <- as.character(args[6])
cores <- as.numeric(args[7])
umap_plot_1_path <- as.character(args[8])
umap_plot_2_path <- as.character(args[9])
umap_tsv_path <- as.character(args[10])

# Load files and format 

con_fasta <- readDNAStringSet(fasta_file)

align_fasta <- readDNAStringSet(align_file)

dist_dt <- fread(dist_file, stringsAsFactors = FALSE)
dist_mat <- as.matrix(dist_dt[, -1, with = FALSE])
colnames(dist_mat) <- colnames(dist_dt[, -1])
rownames(dist_mat) <- as.vector(dist_dt[[1]])
rm(dist_dt)
gc()

meta_dt <- fread(meta_file, stringsAsFactors = FALSE)
meta_df <- as.data.frame(meta_dt)
meta_df_sub <- meta_df[,c("Accession", "Region", "Geo_Location", "Datetime")]
colnames(meta_df_sub) <- c("Accession", "Region", "Country", "Date")
meta_travel_sub <- meta_df[,c("Accession", "Country_Exposure")]
umap_dt <- fread(umap_file, stringsAsFactors = FALSE)
umap_df <- as.data.frame(umap_dt)

print("Step 1 complete - files loaded and reformatted")

# Compute updated alignment

align_new <- align_get(con_fasta, align_fasta)

print("Step 2 complete - alignment updated with new fasta")

# Compute updated distance matrix

dist_mat_up <- dist_get_heur(align_new, con_fasta, dist_mat)

print("Step 3 complete - distance updated with new fasta")

# Process metadata with SIGNAL output

new_accessions <- rownames(dist_mat_up)[(nrow(meta_df_sub)+1):nrow(dist_mat_up)]
meta_new <- data.frame("Accession" = new_accessions, "Region" = paste("Novel", seq(1, length(new_accessions), 1)), "Country" = paste("Novel", seq(1, length(new_accessions), 1)), "Date" = rep((unclass(Sys.Date()) - unclass(as.Date("2019-12-01", format = "%Y-%m-%d"))), length(new_accessions)))
new_meta <- rbind(meta_df_sub, meta_new)

print("Step 4 complete - metadata updated with new fasta")

# Get updated umap data

umap_new <- umap_process_heur(align_new, con_fasta, dist_map_up, new_meta, umap_df)

print(umap_new)

print("Step 5 complete - umap data reprocessed with new fasta")

# Zoom into umap coords and subset

umap_coords_1 <- umap_new[nrow(umap_new), "UMAP_1"]
umap_coords_2 <- umap_new[nrow(umap_new), "UMAP_2"]

coords_1_high <- umap_coords_1 + 5
coords_1_low <- umap_coords_1 - 5
coords_2_high <- umap_coords_2 + 5
coords_2_low <- umap_coords_2 - 5

umap_sub <- umap_new[which((umap_new$UMAP_1 > coords_1_low) & (umap_new$UMAP_1 < coords_1_high) & (umap_new$UMAP_2 > coords_2_low) & (umap_new$UMAP_2 < -15)), ]

print("Step 6 complete - umap data subsetted for proximity to fasta")

# Create and save visualizations

ggplot(data = umap_new, aes(x = UMAP_1, y = UMAP_2)) +
	    theme_few() +
	    geom_jitter(aes(fill = Region), size = 3, position = "jitter", colour = "black", pch = 21, stroke = 0.25) +
	    scale_fill_manual(name = "", values = kev_palette[1:length(unique(umap_new$Region))]) +
	    labs(x = "UMAP 1", y = "UMAP 2") +
	    theme(axis.ticks.x = element_blank()) +
	    theme(axis.ticks.y = element_blank()) +
	    theme(axis.text.y = element_blank()) +
	    theme(axis.text.x = element_blank()) +
	    theme(axis.title.y = element_text(size = 16, face = "bold")) +
	    theme(axis.title.x = element_text(size = 16, face = "bold")) +
	    theme(legend.title = element_text(size = 15, face = "bold")) +
	    theme(legend.text = element_text(size = 14)) 


ggsave(umap_plot_1_path, device = "pdf", width = 12, height = 6)

ggplot(data = umap_sub, aes(x = UMAP_1, y = UMAP_2)) +
       theme_few() +
       geom_jitter(aes(fill = Country), size = 3, position = "jitter", colour = "black", pch = 21, stroke = 0.25) +
       scale_fill_manual(name = "", values = kev_palette[1:length(unique(umap_sub$Country))]) +
       labs(x = "UMAP 1", y = "UMAP 2") +
       theme(axis.ticks.x = element_blank()) +
       theme(axis.ticks.y = element_blank()) +
       theme(axis.text.y = element_blank()) +
       theme(axis.text.x = element_blank()) +
       theme(axis.title.y = element_text(size = 16, face = "bold")) +
       theme(axis.title.x = element_text(size = 16, face = "bold")) +
       theme(legend.title = element_text(size = 15, face = "bold")) +
       theme(legend.text = element_text(size = 14)) 

ggsave(umap_plot_2_path, device = "pdf", width = 12, height = 6)

print("Step 7 complete - umap plots saved")

# Output umap dataframe

umap_sub_seq <- umap_sub[nrow(umap_sub), ]
umap_sub_minus <- umap_sub[-nrow(umap_sub), ]

umap_plus_travel <- merge(umap_sub_minus, meta_travel_sub, by = "Accession")
colnames(umap_plus_travel)[ncol(umap_plus_travel)] <- "Travel_History"
umap_plus_travel$Travel_History <- ifelse((umap_plus_travel$Travel_History == umap_plus_travel$Country), "Not available", umap_plus_travel$Travel_History)

umap_sub_seq$Travel_History <- sample_hist

umap_final_df <- rbind(umap_plus_travel, umap_sub_seq)
umap_final_dt <- as.data.table(umap_final_df)
fwrite(umap_final_dt, umap_tsv_path, quote = FALSE, col.names = TRUE, row.names = FALSE) 

print("Step 8 complete - umap data saved as tsv")

