create_training_testing = function(file_path, p, dbSNP_only = TRUE) {
  cat("- Reading file\n")
  mt = read.table(file_path, 
                  sep = "\t", 
                  header = TRUE, 
                  fill = TRUE, 
                  stringsAsFactors = FALSE)
  # remove the dbSNP only
  if (dbSNP_only) mt = mt[mt$id2 != 0, ]
  
  # divide the samples into training and testing sets
  mt_list = unique(mt$case)
  
  no_testing = FALSE
  if (length(mt_list) == p) no_testing = TRUE
  
  cat("- Training\n")
  set.seed(23)
  training_list = sample(mt_list, p)
  training = mt[ mt$case %in% training_list, ]
  training$num_mut_gene = num_mut_gene(training)
  training$case_mut_gene = case_mut_gene(training)
  training = cleaning(training)
  
  if (! no_testing) {
    cat("- Testing\n")
    testing_list = mt_list[! mt_list %in% training_list]
    testing = mt[ mt$case %in% testing_list, ]
    testing$num_mut_gene = num_mut_gene(testing)
    testing$case_mut_gene = case_mut_gene(testing)
    testing = cleaning(testing)
    
    y = list(training, testing)
  } else {
    y = list(training, FALSE)
  }
  
  return(y)
}