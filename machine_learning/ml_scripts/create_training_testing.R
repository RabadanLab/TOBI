create_training_testing = function(file_path, p, dbSNP_only = TRUE, indel = FALSE) {
  cat("- Reading file\n")
  y = read.table(file_path, 
                  sep = "\t", 
                  header = TRUE, 
                  fill = TRUE, 
                  stringsAsFactors = FALSE)
  # remove the dbSNP only
  if (dbSNP_only) y = y[y$id2 != 0, ]
  
  # divide the samples into training and testing sets
  y_list = unique(y$case)
  
  no_testing = FALSE
  if (length(y_list) == p) no_testing = TRUE
  
  cat("- Training\n")
  set.seed(23)
  training_list = sample(y_list, p)
  training = y[ y$case %in% training_list, ]
  training$num_mut_gene = num_mut_gene(training)
  training$case_mut_gene = case_mut_gene(training)
  training = cleaning(training, indel)
  
  if (! no_testing) {
    cat("- Testing\n")
    testing_list = y_list[! y_list %in% training_list]
    testing = y[ y$case %in% testing_list, ]
    testing$num_mut_gene = num_mut_gene(testing)
    testing$case_mut_gene = case_mut_gene(testing)
    testing = cleaning(testing, indel)
    
    y = list(training, testing)
  } else {
    y = list(training, FALSE)
  }
  
  return(y)
}