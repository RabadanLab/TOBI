cleaning = function(x) {
  # this function gets the mutation list AFTER pre-processing step, anc "cleans" it for machine learning step
  # x is the dataframe
  # the output is the cleaned dataframe
  y = subset(x, select = c(Y, case, chrom, pos, id, ref, alt, gene_name, dp, effect, amino_acid_change, warnings,
                           id2, qual, freq, amino_acid_length, gene_freq, gene_length, cnv, tot_cosm_samp, num_mut_gene, case_mut_gene, mq,
                           effect_impact, cosmic_nsamp, pval_som, dbNSFP_1000Gp1, dbNSFP_CADD_phred,
                           dbNSFP_MutationAssessor_score, dbNSFP_MutationTaster_score,
                           dbNSFP_FATHMM_score, dbNSFP_Polyphen2_HDIV_score, dbNSFP_Polyphen2_HVAR_score,
                           dbNSFP_RadialSVM_score, dbNSFP_SIFT_score))
  
  y$Y[y$Y == 0] = "non_som"
  y$Y[y$Y == 1] = "som"
  y$Y = as.factor(y$Y)
  y$Y = factor(y$Y, levels = c("som", "non_som"))
  
  y$id2 = as.factor(y$id2)
  y$effect_impact = as.factor(y$effect_impact)
  
  y$cnv = as.numeric(y$cnv)
  y$cnv[is.na(y$cnv)] = 0
  
  y$warnings = as.numeric(y$warnings)
  y$qual = as.numeric(y$qual)
  
  y$amino_acid_length = as.numeric(y$amino_acid_length)
  y$amino_acid_length[is.na(y$amino_acid_length)] = -1
  
  y$dbNSFP_MutationAssessor_score = as.numeric(y$dbNSFP_MutationAssessor_score)
  y$dbNSFP_MutationAssessor_score[is.na(y$dbNSFP_MutationAssessor_score)] = -1
  
  y$dbNSFP_MutationTaster_score = as.numeric(y$dbNSFP_MutationTaster_score)
  y$dbNSFP_MutationTaster_score[is.na(y$dbNSFP_MutationTaster_score)] = -1
  
  y$cosmic_nsamp = as.numeric(y$cosmic_nsamp)
  y$cosmic_nsamp[is.na(y$cosmic_nsamp)] = -1
  
  y$dbNSFP_1000Gp1 = as.factor(y$dbNSFP_1000Gp1)
  
  y$dbNSFP_CADD_phred = as.numeric(y$dbNSFP_CADD_phred)
  y$dbNSFP_CADD_phred[is.na(y$dbNSFP_CADD_phred)] = -1
  
  y$dbNSFP_FATHMM_score = as.numeric(y$dbNSFP_FATHMM_score)
  y$dbNSFP_FATHMM_score[is.na(y$dbNSFP_FATHMM_score)] = -1
  
  y$dbNSFP_Polyphen2_HDIV_score = as.numeric(y$dbNSFP_Polyphen2_HDIV_score)
  y$dbNSFP_Polyphen2_HDIV_score[is.na(y$dbNSFP_Polyphen2_HDIV_score)] = -1
  
  y$dbNSFP_Polyphen2_HVAR_score = as.numeric(y$dbNSFP_Polyphen2_HVAR_score)
  y$dbNSFP_Polyphen2_HVAR_score[is.na(y$dbNSFP_Polyphen2_HVAR_score)] = -1
  
  y$dbNSFP_RadialSVM_score = as.numeric(y$dbNSFP_RadialSVM_score)
  y$dbNSFP_RadialSVM_score[is.na(y$dbNSFP_RadialSVM_score)] = -1
  
  y$dbNSFP_SIFT_score = as.numeric(y$dbNSFP_SIFT_score)
  y$dbNSFP_SIFT_score[is.na(y$dbNSFP_SIFT_score)] = -1
  return(y)
}