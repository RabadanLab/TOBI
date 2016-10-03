#!/usr/bin/env Rscript

main = function(filepath, outputpath) {
  cat("file path is:", "\n")
  cat(filepath, "\n")
  
  header = readLines(filepath, n = 1)
  if (grepl("#", header)) {
    header = gsub("#", "", header)
    header = gsub("\t\t", "\tX\tX.1", header)
    header = unlist(strsplit(header, split = "\t"))
    mt = read.table(filepath,
                    sep = "\t", 
                    header = FALSE, 
                    quote = "",
                    stringsAsFactors = FALSE,
                    skip = 1)
    colnames(mt) = header
  } else {
    mt = read.table(filepath,
                    sep = "\t", 
                    header = TRUE,
                    quote = "",
                    stringsAsFactors = FALSE)
  }
  
  # Correcting the indels position to be compatible with Cbio
  mt$pos = mt$pos + as.numeric(nchar(mt$ref) > nchar(mt$alt))
  
  print(dim(mt))
  # Applying techn filter
  mt$qual = as.numeric(mt$qual)
  mt$dp4_3 = as.numeric(mt$dp4_3)
  mt$dp4_4 = as.numeric(mt$dp4_4)
  mt$pv4_1 = as.numeric(mt$pv4_1)
  mt$pv4_3 = as.numeric(mt$pv4_3)
  mt$pv4_4 = as.numeric(mt$pv4_4)
  
  print(dim(mt))
  mt = mt[! (is.na(mt$pv4_1) | 
               is.na(mt$pv4_3) | 
               is.na(mt$pv4_4)), ]
  mt = mt[(mt$qual > 60 | 
             (mt$dp4_3 >= 10 & mt$dp4_4 >= 10)) & 
            mt$pv4_1 > 0.01 & 
            mt$pv4_3 > 0.01 & 
            mt$pv4_4 > 0.01, ]
  
  print(dim(mt))
  # Applying biol filter
  mt = mt[mt$common != "1" &
            mt$g5a == "." & 
            mt$g5 == "." &
	    mt$meganormal_id == "." &
            mt$effect != "." &
            mt$effect != "INTRAGENIC" &
            mt$effect != "EXON" &
            mt$effect != "SPLICE_SITE_REGION", ]
  
  print(dim(mt))
  # To remove the src from old cases
  if ("src" %in% colnames(mt)) {
    mt$src = NULL
  }
  
  print(dim(mt))
  # Removing redundant columns
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

