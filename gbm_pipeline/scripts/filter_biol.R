#!/usr/bin/env Rscript

main = function(filepath, outputpath) {
  mt = read.table(filepath,
                  sep = "\t", 
                  header = TRUE, 
                  stringsAsFactors = FALSE)

  mt = mt[mt$common != "1" &
            mt$g5a == "." & 
            mt$g5 == "." &
            mt$meganormal_id == "." &
            mt$effect != "." &
            mt$effect != "INTRAGENIC" &
            mt$effect != "EXON" &
            mt$effect != "SPLICE_SITE_REGION", ]
  
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

