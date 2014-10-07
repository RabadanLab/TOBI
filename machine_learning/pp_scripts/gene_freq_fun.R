gene_freq_fun = function(y, freq) {
  # This function adds gene_freq and gene_length columns to mutations data frame
  freq$freq = with(freq, count / amino_acid_length)
  freq$freq = sapply(freq$freq, function(x) ifelse(is.nan(x),0,x))
  
  tmp_output = foreach (i=1:dim(y)[1], .combine = rbind) %dopar% {
    row_num = which(freq$gene_name == y$gene_name[i])
    if (length(row_num)) {
      c(freq$freq[row_num], freq$length[row_num])
    } else {
      c(NA, NA)
    }
  }
  y$gene_freq = tmp_output[,1]
  y$gene_length = tmp_output[,2]
  y = y[! is.na(y$gene_freq), ]
  
  y
}