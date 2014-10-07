cleaning = function(x, indel = FALSE) {
  # this function gets the mutation list AFTER pre-processing step, 
  # and "cleans" it for machine learning step
  # x is the dataframe
  # the output is the cleaned dataframe
  
  param = c("Y", "case", "chrom", "pos", "id", "ref", "alt", "gene_name", "dp", "effect",
            "amino_acid_change", "warnings", "id2", "qual", "freq", "amino_acid_length",
            "gene_freq", "gene_length", "tot_cosm_samp", "num_mut_gene", "case_mut_gene",
            "mq", "effect_impact", "cosmic_nsamp", "pval_som", "dbNSFP_1000Gp1",
            "dbNSFP_CADD_phred", "dbNSFP_MutationAssessor_score",
            "dbNSFP_MutationTaster_score", "dbNSFP_FATHMM_score",
            "dbNSFP_Polyphen2_HDIV_score", "dbNSFP_Polyphen2_HVAR_score",
            "dbNSFP_RadialSVM_score", "dbNSFP_SIFT_score")
  # add indel related features, if true
  if (indel) {
    param = c(param, "indel", "is_1", "is_2", "mega")
  }
  y = subset(x, select = param)
  
  y$Y = as.factor(y$Y)
  y$Y = factor(y$Y, levels = c("som", "non_som"))
  
  y$id2 = as.factor(y$id2)
  y$effect_impact = as.factor(y$effect_impact)
  
  y$warnings = as.numeric(y$warnings)
  y$qual = as.numeric(y$qual)
  
  columns = c("amino_acid_length", "cosmic_nsamp")
  
  if (indel) {
    columns = c(columns, "indel", "is_1", "is_2")
    for (column in columns) y[, column] = na_minus1(y[, column])
  } else {
    for (column in columns) y[, column] = na_minus1(y[, column])
  }
  
  y
}