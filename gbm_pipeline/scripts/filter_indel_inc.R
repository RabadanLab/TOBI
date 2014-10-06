#!/usr/bin/env Rscript

main = function(filepath, outputpath) {
  mt = read.table(filepath,
                  sep = "\t", 
                  header = TRUE, 
                  stringsAsFactors = FALSE)
  
  mt$ref_len = sapply(mt$ref, nchar)
  mt$alt_len = sapply(mt$alt, nchar)
  
  mt$del = as.numeric(mt$ref_len > mt$alt_len)
  mt$pos = mt$pos + mt$del

  mt$ref_len = NULL
  mt$alt_len = NULL
  mt$del = NULL
  if (length(grep("^X", colnames(mt)))) {
    mt = mt[, - grep("^X", colnames(mt))]
  }
  
  write.table(mt, 
              outputpath, 
              sep = "\t", 
              quote = FALSE, 
              na = ".", 
              row.names = FALSE)
}

if (!interactive()){
  args = commandArgs(TRUE)
  filepath = args[[1]]
  outputpath = args[[2]]
  main(filepath, outputpath)
}

