cleaning_features = function(mt) {
    mt$effect_impact[mt$effect_impact == "."] = 0
    mt$effect_impact[mt$effect_impact == "LOW"] = 1
    mt$effect_impact[mt$effect_impact == "MODERATE"] = 2
    mt$effect_impact[mt$effect_impact == "MODIFIER"] = 3
    mt$effect_impact[mt$effect_impact == "HIGH"] = 4
    
    mt$gt_1[mt$gt_1 == '0/1'] = 0
    mt$gt_1[mt$gt_1 == '1/1'] = 1
    
    mt$sao[mt$sao == "."] = "NA"
    mt$sao[mt$sao == "0"] = "Unspecified"
    mt$sao[mt$sao == "1"] = "Germline"
    
    mt$warnings[mt$warnings != "."] = 1
    mt$warnings[mt$warnings == "."] = 0
    
    mt$som_counts = as.numeric(as.character(mt$som_counts))
    mt$rpb = as.numeric(as.character(mt$rpb))
    mt$exon_rank = as.numeric(as.character(mt$exon_rank))
    mt$amino_acid_length = as.numeric(as.character(mt$amino_acid_length))
    mt$ssr = as.numeric(as.character(mt$ssr))
    mt$fq = as.numeric(mt$fq)
    
    mt$dbNSFP_CADD_phred = comma_to_num(mt$dbNSFP_CADD_phred)
    mt$dbNSFP_ESP6500_AA_AF = comma_to_num(mt$dbNSFP_ESP6500_AA_AF)
    mt$dbNSFP_phastCons100way_vertebrate = comma_to_num(mt$dbNSFP_phastCons100way_vertebrate)
    mt$dbNSFP_1000Gp1_EUR_AF = comma_to_num(mt$dbNSFP_1000Gp1_EUR_AF)
    mt$dbNSFP_1000Gp1 = as.numeric(! is.na(mt$dbNSFP_1000Gp1_EUR_AF))
    mt$dbNSFP_ESP6500_EA_AF = comma_to_num(mt$dbNSFP_ESP6500_EA_AF)
    mt$dbNSFP_FATHMM_score = comma_to_num(mt$dbNSFP_FATHMM_score)
    mt$dbNSFP_LRT_Omega = comma_to_num(mt$dbNSFP_LRT_Omega)
    mt$dbNSFP_LRT_score = comma_to_num(mt$dbNSFP_LRT_score)
    mt$dbNSFP_LR_score = comma_to_num(mt$dbNSFP_LR_score)
    mt$dbNSFP_MutationAssessor_score = comma_to_num(mt$dbNSFP_MutationAssessor_score)
    mt$dbNSFP_MutationTaster_score = comma_to_num(mt$dbNSFP_MutationTaster_score)
    mt$dbNSFP_Polyphen2_HDIV_score = comma_to_num(mt$dbNSFP_Polyphen2_HDIV_score)
    mt$dbNSFP_Polyphen2_HVAR_score = comma_to_num(mt$dbNSFP_Polyphen2_HVAR_score)
    mt$dbNSFP_RadialSVM_score = comma_to_num(mt$dbNSFP_RadialSVM_score)
    mt$dbNSFP_Reliability_index = comma_to_num(mt$dbNSFP_Reliability_index)
    mt$dbNSFP_SIFT_score = comma_to_num(mt$dbNSFP_SIFT_score)
    mt$dbNSFP_SiPhy_29way_logOdds = comma_to_num(mt$dbNSFP_SiPhy_29way_logOdds)
    mt$dbNSFP_phastCons46way_placental = comma_to_num(mt$dbNSFP_phastCons46way_placental)
    mt$dbNSFP_phastCons46way_primate = comma_to_num(mt$dbNSFP_phastCons46way_primate)
    mt$dbNSFP_phyloP100way_vertebrate = comma_to_num(mt$dbNSFP_phyloP100way_vertebrate)
    mt$dbNSFP_phyloP46way_placental = comma_to_num(mt$dbNSFP_phyloP46way_placental)
    mt$dbNSFP_phyloP46way_primate = comma_to_num(mt$dbNSFP_phyloP46way_primate)
    
    mt
}