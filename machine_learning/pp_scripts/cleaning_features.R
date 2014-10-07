cleaning_features = function(y) {
    y$effect_impact[y$effect_impact == "."] = 0
    y$effect_impact[y$effect_impact == "LOW"] = 1
    y$effect_impact[y$effect_impact == "MODERATE"] = 2
    y$effect_impact[y$effect_impact == "MODIFIER"] = 3
    y$effect_impact[y$effect_impact == "HIGH"] = 4
    
    y$gt_1[y$gt_1 == '0/1'] = 0
    y$gt_1[y$gt_1 == '1/1'] = 1
    
    y$sao[y$sao == "."] = "NA"
    y$sao[y$sao == "0"] = "Unspecified"
    y$sao[y$sao == "1"] = "Germline"
    
    y$warnings[y$warnings != "."] = 1
    y$warnings[y$warnings == "."] = 0
    
    to_num_col = c("som_counts", "exon_rank", "amino_acid_length")
    y[to_num_col] = lapply(y[to_num_col], as.numeric)
    
    dbNSFP = c("dbNSFP_CADD_phred", "dbNSFP_ESP6500_AA_AF", 
               "dbNSFP_phastCons100way_vertebrate", "dbNSFP_1000Gp1_EUR_AF", 
               "dbNSFP_ESP6500_EA_AF", "dbNSFP_FATHMM_score", 
               "dbNSFP_LRT_Omega", "dbNSFP_LRT_score", "dbNSFP_LR_score", 
               "dbNSFP_MutationAssessor_score", "dbNSFP_MutationTaster_score", 
               "dbNSFP_Polyphen2_HDIV_score", "dbNSFP_Polyphen2_HVAR_score", 
               "dbNSFP_RadialSVM_score", "dbNSFP_Reliability_index", 
               "dbNSFP_SIFT_score", "dbNSFP_SiPhy_29way_logOdds", 
               "dbNSFP_phastCons46way_placental", "dbNSFP_phastCons46way_primate", 
               "dbNSFP_phyloP100way_vertebrate", "dbNSFP_phyloP46way_placental", 
               "dbNSFP_phyloP46way_primate", "dbNSFP_GERP_NR", "dbNSFP_GERP_RS")
    
    #y[dbNSFP] = lapply(y[dbNSFP], comma_to_num)
    
    y[dbNSFP] = foreach(i=1:dim(y)[1], .combine = rbind) %dopar% {
      comma_to_num(y[i, dbNSFP])
    }
    
    y$dbNSFP_1000Gp1 = ifelse(y$dbNSFP_1000Gp1_EUR_AF == -1, 0, 1)
    y$dbNSFP_1000Gp1_EUR_AF = NULL
    
    y
}