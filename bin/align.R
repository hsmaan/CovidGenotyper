library(Biostrings)
library(DECIPHER)
library(stringi)
library(stringr)

align_get <- function(fasta, align) {
  
  covid_seq <- readDNAStringSet(fasta)
  covid_align <- align
  dec_align <- AlignProfiles(covid_seq, covid_align, restrict = c(-100, 2, 10))
  fasta_mat <- as.matrix(dec_align)
  fasta_bin <- as.DNAbin(fasta_mat)
  fasta_ungapped <- del.colgapsonly(fasta_bin, threshold = 0.95)
  fasta_string <- fasta_ungapped %>% as.list %>% as.character %>% lapply(., paste0, collapse = "") %>% unlist %>% DNAStringSet
  fasta_final <- subseq(fasta_string, start = 265, end = 29674)
  return(fasta_final)
  
}