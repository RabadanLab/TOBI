create_training_testing = function(file_path, p) {
  cat("- Reading file\n")
  mt = read.table(file_path, 
                  sep = "\t", 
                  header = TRUE, 
                  fill = TRUE, 
                  stringsAsFactors = FALSE)
  # remove the dbSNP only
  mt = mt[mt$id2 != 0, ]
  
  # divide the samples into training and testing sets
  mt_list = unique(mt$case)
  
  cat("- Dividing data into training and testing\n")
  set.seed(23)
  training_list = sample(mt_list, p)
  testing_list = mt_list[! mt_list %in% training_list]
  
  training = mt[ mt$case %in% training_list, ]
  testing = mt[ mt$case %in% testing_list, ]
  
  cat("- Creating num_mut_gene and case_mut_gene\n")
  training$num_mut_gene = num_mut_gene(training)
  training$case_mut_gene = case_mut_gene(training)
  
  testing$num_mut_gene = num_mut_gene(testing)
  testing$case_mut_gene = case_mut_gene(testing)
  
  cat("- Cleaning training and testing\n")
  training = cleaning(training)
  testing = cleaning(testing)
  
  y = list(training, testing)
  return(y)
}