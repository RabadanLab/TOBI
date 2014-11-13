#!/usr/bin/env Rscript

main = function(filepath, outputpath) {
  mt = read.table(filepath,
                  sep = "\t", 
                  header = TRUE, 
                  stringsAsFactors = FALSE)
  
  mt = mt[mt$indel == ".", ]
  
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

