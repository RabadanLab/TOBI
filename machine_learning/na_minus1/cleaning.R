cleaning = function(x, param, indel = FALSE) {
  # this function gets the mutation list AFTER pre-processing step, 
  # and "cleans" it for machine learning step
  # x is the dataframe
  # the output is the cleaned dataframe
  
  # add indel related features, if true
  if (indel) {
    param = c(param, "indel", "is_1", "is_2", "mega")
  }
  # test if parameters match with column names
  x_colnames = colnames(x)
  if (any(! param %in% x_colnames)) {
    not_found = param[! param %in% x_colnames]
    stop(paste(not_found, collapse = " & "), " not found. Check if all 
         parameters (param) are present as columns in mutation df.")
  }
  
  y = subset(x, select = param)
  
  y$Y = as.factor(y$Y)
  y$Y = factor(y$Y, levels = c("som", "non_som"))
  
  y$id2 = as.factor(y$id2)
  y$effect_impact = as.factor(y$effect_impact)
  
  y$warnings = as.numeric(y$warnings)
  if ("qual" %in% colnames(y)) y$qual = as.numeric(y$qual)
  
  columns = c("amino_acid_length", "cosmic_nsamp")
  
  if (indel) {
    columns = c(columns, "indel", "is_1", "is_2")
    for (column in columns) y[, column] = na_to_0(y[, column])
  } else {
    for (column in columns) y[, column] = na_to_0(y[, column])
  }

  dbNSFP = c("dbNSFP_CADD_phred", "dbNSFP_RadialSVM_score", 
               "dbNSFP_MutationAssessor_score", "dbNSFP_MutationTaster_score", 
               "dbNSFP_Polyphen2_HDIV_score", "dbNSFP_Polyphen2_HVAR_score", 
               "dbNSFP_SIFT_score", "dbNSFP_FATHMM_score")
  for (column in dbNSFP) y[, column] = na_to_mean(y[, column]) 
  y
}
