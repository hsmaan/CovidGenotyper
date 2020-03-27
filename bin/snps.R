library(reshape2)
library(dplyr)

snps_get <- function(alignment, metadata) {
  
  align <- alignment
  pos <- c(1059,1190, 3037, 17858, 18060, 23403, 25563, 27046)
  acc <- metadata[,1]
  meta_var <- metadata[,2]
  align_df <- as.data.frame(as.matrix(align))
  meta_order <- match(rownames(align_df), acc)
  align_df <- align_df[order(meta_order),]
  align_df$meta <- meta_var
  
  freq_pct <- function(d_col) {
    freq_table <- table(d_col)
    pct_table <- as.table(sapply(freq_table, function(x) x/sum(freq_table)))
    return(pct_table)
  }
  
  table_get <- function(df, position) {
    align_pos <- df[,c(position, ncol(df))]
    align_grouped <- group_by(.data = align_pos, meta) 
    align_tables <- as.data.frame(do(.data = align_grouped, data.frame(val = freq_pct(.[,1]))))
    align_final <- data.frame("Position" = rep(position, length(align_tables[,1])), "Meta" = align_tables[,1], "Allele" = align_tables[,2], "Freq" = align_tables[,3])
    return(align_final)
  }

  all_tables <- lapply(pos, function(x) table_get(align_df, x))
  table_concat <- base::Reduce(rbind, all_tables)
  return(table_concat)
  
}


