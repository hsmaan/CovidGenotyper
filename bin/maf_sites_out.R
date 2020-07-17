library(data.table)
library(tidyverse)
library(Biostrings)
library(ape)
library(GenomicRanges)
library(stringr)

# Load data

setwd("../data")

loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

file_list <- list.files()

alignment <- readDNAStringSet(grep("dec_aligned_fasta_filtered*", list.files(), value = TRUE))
gff <- read.gff(grep("ncov_NC_045512_Genes.GFF3", list.files(), value = TRUE))
var_freq <- read.table(grep("*.frq", list.files(), value = TRUE), stringsAsFactors = FALSE, fill = TRUE, col.names = paste0("V", seq(0, 15)))
vcf <- fread(grep("*.vcf", list.files(), value = TRUE))
vcf <- vcf[-1,]

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

var_freq_filtered <- var_freq_filtered[which((var_freq_filtered$A1 %in% c("A", "T", "C", "G")) & (var_freq_filtered$A2 %in% c("A", "T", "C", "G"))), ]

# Format vcf metadata and add 

vcf_sub <- data.frame("Pos" = vcf$POS, "Meta" = str_split_fixed(vcf$INFO, stringr::fixed("|"), 4)[,2])

var_freq_filtered <- merge(var_freq_filtered, vcf_sub, by = "Pos")

# Covert to GRanges object

var_ranges <- GRanges(seqnames = rep(1, length(var_freq_filtered$Pos)), ranges = IRanges(start = var_freq_filtered$Pos, end = var_freq_filtered$Pos))

var_ranges$effect <- var_freq_filtered$Meta

# Overlap var freq ranges with annotations and remerge

var_overlap <- as.data.frame(mergeByOverlaps(var_ranges, gff_ranges))

var_overlap_sub <- var_overlap[, c("var_ranges.start", "gene", "effect")]
colnames(var_overlap_sub) <- c("Pos", "Gene", "Effect")
var_overlap_sub$Effect <- gsub("_", " ", var_overlap_sub$Effect)

var_freq_overlap <- merge(var_freq_filtered, var_overlap_sub)

var_freq_overlap <- var_freq_overlap[,c("Pos", "A1", "AF1", "A2", "AF2", "Gene", "Effect")]

var_freq_overlap <- var_freq_overlap[order(var_freq_overlap$AF2, decreasing = TRUE),]

# Output data

file.remove(grep("var_freq_*", file_list, value = TRUE))
save(var_freq_overlap, file = paste("var_freq_", Sys.Date(), ".RData", sep = ""))  
save(var_freq_overlap, file = paste("signal_var_freq_", Sys.Date(). ".RData", sep = "")) 
