library(data.table)
library(tidyverse)
library(Biostrings)
library(ape)
library(GenomicRanges)

# Load data

setwd("../data")

loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

file_list <- list.files()

alignment <- readDNAStringSet(grep("dec_aligned_fasta_filtered*", list.files(), value = TRUE))
fasta_ref <- readDNAStringSet(grep("ncov_ref_NC_045512.fasta", list.files(), value = TRUE))
gff <- read.gff(grep("ncov_NC_045512_Genes.GFF3", list.files(), value = TRUE))
var_freq <- read.table(grep("*.frq", list.files(), value = TRUE), stringsAsFactors = FALSE, fill = TRUE, col.names = paste0("V", seq(0, 10)))

# Format and Convert GFF to GRanges object 

gff_gene_sub <- gff[which(gff$type == "gene"),]
gff_gene_names <- str_split_fixed(gff_gene_sub$attributes, "=", 6)
gff_gene_names <- str_split_fixed(gff_gene_names[, 5], ";", 2)[,1]
gff_gene_sub$gene_name <- gff_gene_names

gff_ranges <- GRanges(seqnames = rep(1, length(gff_gene_sub$seqid)), ranges = IRanges(start = gff_gene_sub$start, end = gff_gene_sub$end))
gff_ranges$gene <- gff_gene_sub$gene_name

# Reformat columns for var_freq

var_freq <- as.data.frame(var_freq[-1,])

var_freq_pos <- as.numeric(var_freq$V1)
var_freq_mj_a1 <- str_split_fixed(var_freq$V4, ":", 2)[,1]
var_freq_mj_af1 <- as.numeric(str_split_fixed(var_freq$V4, ":", 2)[,2])
var_freq_mj_a2 <- str_split_fixed(var_freq$V5, ":", 2)[,1]
var_freq_mi_af2 <- as.numeric(str_split_fixed(var_freq$V5, ":", 2)[,2])

var_freq_filtered <- data.frame("Pos" = var_freq_pos, "A1" = var_freq_mj_a1, "AF1" = var_freq_mj_af1, "A2" = var_freq_mj_a2, "AF2" = var_freq_mi_af2)

# Covert to GRanges object

var_ranges <- GRanges(seqnames = rep(1, length(var_freq_filtered$Pos)), ranges = IRanges(start = var_freq_filtered$Pos, end = var_freq_filtered$Pos))

var_freq_sub <- var_freq_filtered[which((var_freq_filtered$AF2 >= 0.005) & (var_freq_filtered$AF1 >= 0.005)), ]

var_freq_sub_9 <- (var_freq_sub[order(var_freq_sub$AF2, decreasing = TRUE),])[1:9,]



  