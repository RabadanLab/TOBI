cleaning_features = function(df) {
  df$effect_impact[df$effect_impact == "."] <- 0
  df$effect_impact[df$effect_impact == "LOW"] <- 1
  df$effect_impact[df$effect_impact == "MODERATE"] <- 2
  df$effect_impact[df$effect_impact == "MODIFIER"] <- 3
  df$effect_impact[df$effect_impact == "HIGH"] <- 4
  
  if ("gt_1" %in% colnames(df)) {
    df$gt_1[df$gt_1 == '0/1'] <- 0
    df$gt_1[df$gt_1 == '1/1'] <- 1
  }
  
  df$sao[df$sao == "."] <- "NA"
  df$sao[df$sao == "0"] <- "Unspecified"
  df$sao[df$sao == "1"] <- "Germline"
  
  df$warnings[df$warnings != "."] <- 1
  df$warnings[df$warnings == "."] <- 0
  
  to_num_col <- c("som_counts", "exon_rank", "amino_acid_length")
  df[to_num_col] <- lapply(df[to_num_col], as.numeric)
  
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
  
  #df[dbNSFP] = lapply(df[dbNSFP], comma_to_num)
  
  df[dbNSFP] = foreach(i=1:nrow(df), .combine = rbind) %dopar% {
    comma_to_num(df[i, dbNSFP])
  }
  
  df$dbNSFP_1000Gp1 = ifelse(df$dbNSFP_1000Gp1_EUR_AF == -1, 0, 1)
  df$dbNSFP_1000Gp1_EUR_AF = NULL
  
  return(df)
}

comma_to_num = function(input) {
  # This function cleans a vector with multiple values in random cells,
  # and replace the NAs with -1
  output = sapply(input, function(input) {
    max(as.numeric(unlist(strsplit(as.character(input), ","))))
  })
  output = unname(output)
  output = as.numeric(output)
  output[is.na(output)] = -1
  return(unname(output))
}

gene_freq_fun = function(df, freq) {
  # This function adds gene_freq and gene_length columns to mutations data frame
  freq$freq = with(freq, count / amino_acid_length)
  freq$freq = sapply(freq$freq, function(x) ifelse(is.nan(x),0,x))
  
  tmp_output = foreach (i=1:nrow(df), .combine = rbind) %dopar% {
    row_num = which(freq$gene_name == df$gene_name[i])
    if (length(row_num)) {
      c(freq$freq[row_num], freq$length[row_num])
    } else {
      c(NA, NA)
    }
  }
  df$gene_freq = tmp_output[,1]
  df$gene_length = tmp_output[,2]
  df = df[! is.na(df$gene_freq), ]
  
  return(df)
}

indel_filter = function(df, mega) {
  # flag indels close to meganormals
  l = dim(df)[1]
  df$mega = foreach (i = 1:l, .combine = c) %dopar% {
    mega_fun(df$chrom[i], df$pos[i], mega)
  }
  
  #   table(df$mega[df$indel == 1])
  #   table(df$mega[df$indel == 1 & df$Y == "som"])
  
  df$ref_len = nchar(df$ref)
  df$alt_len = nchar(df$alt)
  
  df$is_1 = as.numeric(df$is_1)
  df$is_2 = as.numeric(df$is_2)
  
  # only deletions
  df = df[ df$indel == "." | (df$ref_len > df$alt_len), ]
  
  # only Frame_shift
  df = df[ df$indel == "." | (df$effect == "FRAME_SHIFT"), ]
  
  # indels longer than 10 are excluded
  df = df[ df$ref_len <= 10 & df$alt_len <= 10, ]
  
  # remvoing indels close to meganormal
  df = df[ df$indel == "." | df$mega == 0, ]
  
  # removing indels with is_1 < 100
  df = df[ df$indel == "." | df$is_1 < 100, ]
  
  #   table(df$Y[df$indel == "1"])
  
  df$ref_len = NULL
  df$alt_len = NULL
  
  return(df)
}

mega_fun = function(chr, pos, mega) {
  # this function determines if a mutation is within the intervals of meganormal
  # mutations, mega
  if (mega$start[mega$chr == chr][1] <= pos) {
    n = which(mega$start[mega$chr == chr] > pos)[1]
    ifelse(mega[n-1, 3] > pos, 1, 0)
  } else {
    0
  }
}

id2 = function(df) {
  # .
  df$in_none = as.numeric(! grepl("rs", df$id) & ! grepl("COS", df$id))*1   
  # only dbSNP
  df$in_dbsnp = as.numeric(grepl("rs", df$id) & ! grepl("COS", df$id))*0
  # only Cosmic
  df$in_cosmic = as.numeric(grepl("COS", df$id) & ! grepl("rs", df$id))*3
  # both Cosmic and dbSNP
  df$in_both = as.numeric(grepl("rs", df$id) & grepl("COS", df$id))*2       
  
  return(df$in_none + df$in_dbsnp + df$in_cosmic + df$in_both)
}

my_read_table = function(x) {
  read.table(x, 
             #sep = "\t", 
             header = TRUE,
             quote = "",
             stringsAsFactors = FALSE)
}

my_write_table = function(variable, output_file) {
  write.table(variable, file = output_file,
              quote = FALSE, 
              sep = "\t", 
              na = ".", 
              row.names = FALSE, 
              col.names = TRUE)
  cat("- File written.\n")
}

removing_features = function(df) {
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
  
  df = df[,!(colnames(df) %in% drops)]
  
  df = df[, -grep("*rankscore*", colnames(df))]
  df = df[, -grep("*pred*", colnames(df))]
  df = df[, -grep("*cln*", colnames(df))]
  
  return(df)
}

pfam_cosmic_fun = function(df, pfam_cosmic) {
  # This function adds the total number of cosmic mutations in a given domain
  # to the mutations that happen within that domain. The rest will have value of 0
  library(doMC)
  registerDoMC(cores = 4)
  
  df$pfam_cosmic = foreach (i = 1:dim(df)[1], .combine = "c") %dopar% {
    ind = (pfam_cosmic$chromStart[ pfam_cosmic$chrom == df$chrom[i]] <= df$pos[i]) &
      (pfam_cosmic$chromEnd[ pfam_cosmic$chrom == df$chrom[i]] >= df$pos[i])
    if (any(ind)) {
      max(pfam_cosmic$pval[ind])
    } else {
      0
    }
  }
  return(df)
}

somatic_annot = function(df, somatic) {
  # somatic_file is a list of somatic mutations with "case_id", "chr", and
  # "start_position"
  
  y_list = unique(df$case)
  som_list = unique(somatic$case_id)
  
  # Choose the corresponding cases from somatic
  somatic = subset(somatic, somatic$case_id %in% y_list)
  
  somatic$ID = paste(somatic$case_id, somatic$chr, 
                     somatic$start_position, somatic$variant_allele, sep = ":")
  df$ID = paste(df$case, df$chrom, df$pos, df$alt, sep = ":")
  
  # Finding the somatic mutations using CBio data
  somatic$Not_Present = 1
  df$Y = "non_som"
  
  for (id in somatic$ID) {
    a = which(df$ID == id)
    if (length(a)) {
      df$Y[a] = "som"
      somatic$Not_Present[somatic$ID == id] = 0
    }
  }
  
  return(list(df, somatic))
}

tot_cosm_fun = function(df, cosm) {
  l = dim(df)[1]
  df$tot_cosm_samp = foreach (i = 1:l, .combine = c) %dopar% {
    row_num = which(cosm$gene == df$gene_name[i])
    if (length(row_num) > 0) {
      cosm$tot_nsample[row_num]
    } else {
      0
    }
  }
  return(df)
}