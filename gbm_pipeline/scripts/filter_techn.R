#!/usr/bin/env Rscript

main = function(filepath, outputpath) {
  mt = read.table(filepath,
                  sep = "\t", 
                  header = TRUE, 
                  stringsAsFactors = FALSE)
  
  mt$qual = as.numeric(mt$qual)
  mt$dp4_3 = as.numeric(mt$dp4_3)
  mt$dp4_4 = as.numeric(mt$dp4_4)
  mt$pv4_1 = as.numeric(mt$pv4_1)
  mt$pv4_3 = as.numeric(mt$pv4_3)
  mt$pv4_4 = as.numeric(mt$pv4_4)
  
  mt = mt[! (is.na(mt$pv4_1) | 
               is.na(mt$pv4_3) | 
               is.na(mt$pv4_4)), ]
  mt = mt[(mt$qual > 60 | 
             (mt$dp4_3 >= 10 & mt$dp4_4 >= 10)) & 
            mt$pv4_1 > 0.01 & 
            mt$pv4_3 > 0.01 & 
            mt$pv4_4 > 0.01, ]
  
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

