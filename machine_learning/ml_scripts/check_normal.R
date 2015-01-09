check_normal = function(file_path, modboost, param) {
  # file_path is the normal data after pre-processing
  cat("- Reading in the normal data\n")
  mt = read.table(file_path, 
                  sep = "\t", 
                  header = TRUE, 
                  fill = TRUE, 
                  stringsAsFactors = FALSE)
  # remove the dbSNP only
  mt = mt[mt$id2 != 0, ]
  
  cat("- Creating num_mut_gene and case_mut_gene\n")
  mt$num_mut_gene = num_mut_gene(mt)
  mt$case_mut_gene = case_mut_gene(mt)
  
  cat("- Cleaning the data\n")
  mt = cleaning(mt, param)
  
  cat("- Confusion table\n")
  print(table(predict(modboost, mt[,-c(1:11)]), mt[,1], dnn=list('predicted','actual')))
}