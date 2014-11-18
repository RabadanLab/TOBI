removing_features = function(y) {
  # This function removes the redundant columns
  drops = c("aa", "ac", "af1", "an", "asp", "ass", "caf", "cda", "cds", 
            "cfl", "cgt", "clnacc", "clndsdb", "clndsdbid", "clnhgvs", 
            "clnsrc", "clnsrcid", "clr", "common", "dbNSFP_1000Gp1_AC", 
            "dbNSFP_1000Gp1_AF", "dbNSFP_1000Gp1_AFR_AC", 
            "dbNSFP_1000Gp1_AFR_AF", "dbNSFP_1000Gp1_AMR_AC", 
            "dbNSFP_1000Gp1_AMR_AF", "dbNSFP_1000Gp1_ASN_AC", 
            "dbNSFP_1000Gp1_ASN_AF", "dbNSFP_1000Gp1_EUR_AC", 
            "dbNSFP_1000Gp1_EUR_AF", "dbNSFP_CADD_raw", 
            "dbNSFP_CADD_raw_rankscore", "dbNSFP_SiPhy_29way_pi", 
            "dbNSFP_UniSNP_ids", "dbSNPBuildID", "dss", "filter", "fq", 
            "g3", "g5", "g5a", "gene", "geneinfo", "genotype_number", "hwe", 
            "ID", "indel.1", "kgvalidated", "lof", "mdv", "meganormal_id", 
            "yp", "mut", "nmd", "nmutperid", "noc", "nov", "nsf", "nsn", 
            "om", "oth", "pc2", "pchi2", "pl_1", "pmc", "pr", "qbd", 
            "qchi2", "r3", "rpq", "rs", "rspos", "ssr", "strand", "tpa", 
            "transcript_id", "u3", "u5", "ugt", "vp", "wtd", "X", "X.1")
  
  y = y[,!(colnames(y) %in% drops)]
   
  y = y[, -grep("*rankscore*", colnames(y))]
  y = y[, -grep("*pred*", colnames(y))]
  y = y[, -grep("*cln*", colnames(y))]
  
  y
}