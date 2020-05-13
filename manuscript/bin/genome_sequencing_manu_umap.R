library(Biostrings)
library(DECIPHER)
library(stringr)
library(ape)
library(reshape2)
library(dplyr)
library(uwot)
library(ggplot2)
library(ggthemes)
library(Cairo)
library(data.table)

# Source all functions

setwd("../../") # main CGT dir

source("R/global.R")

# Load palettes

kev_palette <- c(
  "dodgerblue2", "#E31A1C",
  "green4",
  "#6A3D9A", 
  "#FF7F00", 
  "black", "gold1",
  "skyblue2", "#FB9A99", 
  "palegreen2",
  "#CAB2D6", 
  "#FDBF6F", 
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

qual_palettes = brewer.pal.info[brewer.pal.info$category == "qual", ]
qual_vector = unlist(mapply(brewer.pal, qual_palettes$maxcolors, rownames(qual_palettes)))

# Load aligned GISAID profile and metadata

data_files <- list.files("data")

setwd("data")

loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

pre_aligned <- loadRData(grep("dec_aligned_filtered", data_files, value = TRUE))
pre_meta <- loadRData(grep("covid_filtered_meta", data_files, value = TRUE))
meta <- pre_meta[,c("Accession", "Region", "Geo_Location", "Datetime")]

# Subset pre_aligned data for pre mar 31 (datetime = 121) data

meta_mar31 <- pre_meta[which(meta$Datetime <= 121), ]

pre_aligned_sub <- pre_aligned[c(meta_mar31$Accession)]

# Read fasta for Iran isolate

iran_seq <- readDNAStringSet("Iran1-LN-unicycler-round2.fasta")

setwd("..")

# Concatenate contigs

iran_seq_concat <- DNAStringSet(c(iran_seq[[2]], iran_seq[[1]]))

# Profile alginmeent to GISAID

align_get <- function(stringset, align) { # slight tweak of app function
  
  covid_seq <- stringset
  covid_align <- align
  covid_seq <- RemoveGaps(covid_seq, removeGaps = "all", processors = NULL)
  fasta_final <- AlignProfiles(covid_align, covid_seq)
  return(fasta_final)
  
}

iran_align <- align_get(iran_seq_concat, pre_aligned_sub)

names(iran_align)[length(iran_align)] <- "Novel_Iran1_LN"

writeXStringSet(iran_align, file = "manuscript/data/iran_align_full.fasta")

# Test 

mask_sites <- c(187, 1059, 2094, 3037, 3130, 6990, 8022, 10323, 10741, 11074, 13408, 14786, 19684, 20148, 21137, 24034, 24378, 25563, 26144, 26461, 26681, 28077, 28826, 28854, 29700, 4050, 13402, 11083, 15324, 21575)

mask_sites <- unique(mask_sites)

# Trim off ends and get seq distances

dist_get <- function(align, masked) { # another tweak of an app function
  
  align_mat <- as.matrix(align)
  align_mat_sub <- align_mat[, -masked]
  align_mat_bin <- as.DNAbin(align_mat_sub)
  align_masked <- align_mat_bin %>% as.list %>% as.character %>% lapply(., paste0, collapse = "") %>% unlist %>% DNAStringSet
  align_trim <- subseq(align_masked, start = 316, end = 29674) # Adjusted to 316 start to accomodate missing UTR, considering doing this for all GISAID in future
  dec_dist <- dist.dna(as.DNAbin(align_trim), model = "K80", as.matrix = TRUE, pairwise.deletion = FALSE)
  colnames(dec_dist) <- (str_split_fixed(colnames(dec_dist), fixed("."), 2)[,1])
  rownames(dec_dist) <- (str_split_fixed(rownames(dec_dist), fixed("."), 2)[,1])
  gc()
  return(dec_dist)
  
}

iran_dist <- dist_get(iran_align, mask_sites)

save(iran_dist, file = "manuscript/data/iran_dist_aligned.RData")

# Adjust metadata 

new_meta <- data.frame("Accession" = c("Novel_Iran1_LN"), "Region" = c("Novel_Iran1_LN"), "Geo_Location" = c("Novel_Iran1_LN"), "Datetime" = (unclass(Sys.Date()) - unclass(as.Date("2019-12-01", format = "%Y-%m-%d"))))

meta_updated <- rbind(meta, new_meta)

# Get and save umap data

iran_umap <- umap_process(iran_dist, meta_updated)

fwrite(as.data.table(iran_umap), file = "manuscript/data/iran_umap_full.csv", sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Get plots

umap_plotter <- function(umap_df) { # Function on hand to adjust 
  
  p1 <- ggplot(data = umap_df, aes(x = UMAP_1, y = UMAP_2)) +
    theme_few() +
    geom_jitter(aes(fill = Region), size = 3, position = "jitter", colour = "black", pch = 21, stroke = 0.25) +
    scale_fill_manual(name = "", values = kev_palette[1:length(unique(umap_df$Region))]) +
    geom_point(data = umap_df[grep("Novel", umap_df$Region), ], pch = 21, fill = NA, size = 4, colour = "firebrick1", stroke = 4) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    theme(axis.text.y = element_blank()) +
    theme(axis.text.x = element_blank()) +
    theme(axis.title.y = element_text(size = 16, face = "bold")) +
    theme(axis.title.x = element_text(size = 16, face = "bold")) +
    theme(legend.title = element_text(size = 15, face = "bold")) +
    theme(legend.text = element_text(size = 14)) +
    theme(aspect.ratio = 0.6)
  
  p2 <- ggplot(data = umap_df, aes(x = UMAP_1, y = UMAP_2)) +
    theme_few() +
    geom_jitter(aes(fill = Geo_Location), size = 3, position = "jitter", colour = "black", pch = 21, stroke = 0.25) +
    scale_fill_manual(name = "", values = qual_vector[1:length(unique(umap_df$Geo_Location))]) +
    geom_point(data = umap_df[grep("Novel", umap_df$Geo_Location), ], pch = 21, fill = NA, size = 4, colour = "firebrick1", stroke = 4) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    theme(axis.text.y = element_blank()) +
    theme(axis.text.x = element_blank()) +
    theme(axis.title.y = element_text(size = 16, face = "bold")) +
    theme(axis.title.x = element_text(size = 16, face = "bold")) +
    theme(legend.title = element_text(size = 15, face = "bold")) +
    theme(legend.text = element_text(size = 14)) +
    theme(aspect.ratio = 0.6)
  
  p3 <- ggplot(data = umap_df, aes(x = UMAP_1, y = UMAP_2)) +
    theme_few() +
    geom_jitter(aes(fill = Datetime), size = 3, position = "jitter", colour = "black", pch = 21, stroke = 0.25) +
    scale_fill_gradientn(name = "", colours = c("dodgerblue2", "white", "firebrick1"), breaks = c(min(umap_df$Datetime), median(umap_df$Datetime), max(umap_df$Datetime)), labels = c("Early", "Mid" , "Late")) +
    geom_point(data = umap_df[grep("Novel", umap_df$Geo_Location), ], pch = 21, fill = NA, size = 4, colour = "firebrick1", stroke = 4) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    theme(axis.text.y = element_blank()) +
    theme(axis.text.x = element_blank()) +
    theme(axis.title.y = element_text(size = 16, face = "bold")) +
    theme(axis.title.x = element_text(size = 16, face = "bold")) +
    theme(legend.title = element_text(size = 15, face = "bold")) +
    theme(legend.text = element_text(size = 14)) +
    theme(aspect.ratio = 0.6)
  
  plot_list <- list(p1, p2, p3)
  return(plot_list)
  
}

iran_plots <- umap_plotter(iran_umap)

# View by region

iran_plots[1]

# Sequence falls within a smaller cluster - get umap coords 

ggplot(data = iran_umap, aes(x = UMAP_1, y = UMAP_2)) +
  theme_few() +
  geom_jitter(aes(fill = Region), size = 3, position = "jitter", colour = "black", pch = 21, stroke = 0.25) +
  scale_fill_manual(name = "", values = kev_palette[1:length(unique(iran_umap$Region))]) +
  geom_point(data = iran_umap[grep("Novel", iran_umap$Region), ], pch = 21, fill = NA, size = 4, colour = "firebrick1", stroke = 4) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme(axis.title.y = element_text(size = 16, face = "bold")) +
  theme(axis.title.x = element_text(size = 16, face = "bold")) +
  theme(legend.title = element_text(size = 15, face = "bold")) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(legend.text = element_text(size = 14)) +
  theme(aspect.ratio = 0.6)

ggsave("manuscript/data/figures/Iran1_LN_Raw_UMAP.pdf", device = cairo_pdf, width = 12, height = 6)

# Subset by coords

iran_umap$Region <- as.factor(iran_umap$Region)

iran_umap_sub <- iran_umap[(iran_umap$UMAP_1) > -110 & (iran_umap$UMAP_1 < -60) & (iran_umap$UMAP_2 > -30) & (iran_umap$UMAP_2 < 5), ]

ggplot(data = iran_umap_sub, aes(x = UMAP_1, y = UMAP_2)) +
  theme_few() +
  geom_jitter(aes(fill = Region), size = 3, position = "jitter", colour = "black", pch = 21, stroke = 0.25) +
  scale_fill_manual(name = "", values = kev_palette[1:length(unique(iran_umap$Region))], drop = FALSE) +
  geom_point(data = iran_umap_sub[grep("Novel", iran_umap_sub$Region), ], pch = 21, fill = NA, size = 4, colour = "firebrick1", stroke = 4) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme(axis.title.y = element_text(size = 16, face = "bold")) +
  theme(axis.title.x = element_text(size = 16, face = "bold")) +
  theme(legend.title = element_text(size = 15, face = "bold")) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(legend.text = element_text(size = 14)) +
  theme(aspect.ratio = 0.6)

ggsave("manuscript/data/figures/Iran1_LN_UMAP_Reg_Zoom.pdf", device = cairo_pdf, width = 12, height = 6)

ggplot(data = iran_umap_sub, aes(x = UMAP_1, y = UMAP_2)) +
  theme_few() +
  geom_jitter(aes(fill = Geo_Location), size = 3, position = "jitter", colour = "black", pch = 21, stroke = 0.25) +
  scale_fill_manual(name = "", values = qual_vector[1:length(unique(iran_umap_sub$Geo_Location))]) +
  geom_point(data = iran_umap_sub[grep("Novel", iran_umap_sub$Geo_Location), ], pch = 21, fill = NA, size = 4, colour = "firebrick1", stroke = 4) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme(axis.title.y = element_text(size = 16, face = "bold")) +
  theme(axis.title.x = element_text(size = 16, face = "bold")) +
  theme(legend.title = element_text(size = 15, face = "bold")) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(legend.text = element_text(size = 14)) +
  theme(aspect.ratio = 0.6)

ggsave("manuscript/data/figures/Iran1_LN_UMAP_Country_Zoom.pdf", device = cairo_pdf, width = 12, height = 6)

# Plots no circle 

ggplot(data = iran_umap, aes(x = UMAP_1, y = UMAP_2)) +
  theme_few() +
  geom_jitter(aes(fill = Region), size = 3, position = "jitter", colour = "black", pch = 21, stroke = 0.25) +
  scale_fill_manual(name = "", values = kev_palette[1:length(unique(iran_umap$Region))]) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme(axis.title.y = element_text(size = 16, face = "bold")) +
  theme(axis.title.x = element_text(size = 16, face = "bold")) +
  theme(legend.title = element_text(size = 15, face = "bold")) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(legend.text = element_text(size = 14)) +
  theme(aspect.ratio = 0.6)

ggsave("manuscript/data/figures/Iran1_LN_No_Circle_Raw_UMAP.pdf", device = cairo_pdf, width = 12, height = 6)

ggplot(data = iran_umap_sub, aes(x = UMAP_1, y = UMAP_2)) +
  theme_few() +
  geom_jitter(aes(fill = Region), size = 3, position = "jitter", colour = "black", pch = 21, stroke = 0.25) +
  scale_fill_manual(name = "", values = kev_palette[1:length(unique(iran_umap$Region))], drop = FALSE) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme(axis.title.y = element_text(size = 16, face = "bold")) +
  theme(axis.title.x = element_text(size = 16, face = "bold")) +
  theme(legend.title = element_text(size = 15, face = "bold")) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(legend.text = element_text(size = 14)) +
  theme(aspect.ratio = 0.6)

ggsave("manuscript/data/figures/Iran1_LN_No_Circle_UMAP_Reg_Zoom.pdf", device = cairo_pdf, width = 12, height = 6)

ggplot(data = iran_umap_sub, aes(x = UMAP_1, y = UMAP_2)) +
  theme_few() +
  geom_jitter(aes(fill = Geo_Location), size = 3, position = "jitter", colour = "black", pch = 21, stroke = 0.25) +
  scale_fill_manual(name = "", values = qual_vector[1:length(unique(iran_umap_sub$Geo_Location))]) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme(axis.title.y = element_text(size = 16, face = "bold")) +
  theme(axis.title.x = element_text(size = 16, face = "bold")) +
  theme(legend.title = element_text(size = 15, face = "bold")) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(legend.text = element_text(size = 14)) +
  theme(aspect.ratio = 0.6)

ggsave("manuscript/data/figures/Iran1_LN_No_Circle_UMAP_Country_Zoom.pdf", device = cairo_pdf, width = 12, height = 6)

# Output dataframe and coordinates for iran cluster 

fwrite(as.data.table(iran_umap_sub), file = "manuscript/data/iran_umap_sub.csv", sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Get cluster for iran isolate and output 

iran_umap_cluster <- iran_umap_sub[(iran_umap_sub$UMAP_1) > -87 & (iran_umap_sub$UMAP_1 < -80) & (iran_umap_sub$UMAP_2 > -20) & (iran_umap_sub$UMAP_2 < -12), ]

ggplot(data = iran_umap_cluster, aes(x = UMAP_1, y = UMAP_2)) +
  theme_few() +
  geom_jitter(aes(fill = Region), size = 3, position = "jitter", colour = "black", pch = 21, stroke = 0.25) +
  scale_fill_manual(name = "", values = kev_palette[1:length(unique(iran_umap$Region))]) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  coord_cartesian(xlim = c(-100, -70), ylim = c(-30, 0)) +
  theme(axis.title.y = element_text(size = 16, face = "bold")) +
  theme(axis.title.x = element_text(size = 16, face = "bold")) +
  theme(legend.title = element_text(size = 15, face = "bold")) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(legend.text = element_text(size = 14)) +
  theme(aspect.ratio = 0.6)

fwrite(as.data.table(iran_umap_cluster), file = "manuscript/data/iran_umap_cluster.csv", sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Different palette for iran sub 

ggplot(data = iran_umap_sub, aes(x = UMAP_1, y = UMAP_2)) +
  theme_few() +
  geom_jitter(aes(fill = Geo_Location), size = 3, position = "jitter", colour = "black", pch = 21, stroke = 0.25) +
  scale_fill_manual(name = "", values = kev_palette[1:length(unique(iran_umap_sub$Geo_Location))]) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme(axis.title.y = element_text(size = 16, face = "bold")) +
  theme(axis.title.x = element_text(size = 16, face = "bold")) +
  theme(legend.title = element_text(size = 15, face = "bold")) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(legend.text = element_text(size = 14)) +
  theme(aspect.ratio = 0.6)

ggsave("manuscript/data/figures/Iran1_LN_No_Circle_KevPalette_UMAP_Country_Zoom.pdf", device = cairo_pdf, width = 12, height = 6)
