cleaning_features = function(input) {
  input$effect_impact[input$effect_impact == "."] <- 0
  input$effect_impact[input$effect_impact == "LOW"] <- 1
  input$effect_impact[input$effect_impact == "MODERATE"] <- 2
  input$effect_impact[input$effect_impact == "MODIFIER"] <- 3
  input$effect_impact[input$effect_impact == "HIGH"] <- 4
  
  if ("gt_1" %in% colnames(input)) {
    input$gt_1[input$gt_1 == '0/1'] <- 0
    input$gt_1[input$gt_1 == '1/1'] <- 1
  }
  
  input$sao[input$sao == "."] <- "NA"
  input$sao[input$sao == "0"] <- "Unspecified"
  input$sao[input$sao == "1"] <- "Germline"
  
  input$warnings[input$warnings != "."] <- 1
  input$warnings[input$warnings == "."] <- 0
  
  to_num_col <- c("som_counts", "exon_rank", "amino_acid_length")
  input[to_num_col] <- lapply(input[to_num_col], as.numeric)
  
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
  
  #input[dbNSFP] = lapply(input[dbNSFP], comma_to_num)
  
  input[dbNSFP] = foreach(i=1:nrow(input), .combine = rbind) %dopar% {
    comma_to_num(input[i, dbNSFP])
  }
  
  input$dbNSFP_1000Gp1 = ifelse(input$dbNSFP_1000Gp1_EUR_AF == -1, 0, 1)
  input$dbNSFP_1000Gp1_EUR_AF = NULL
  
  return(input)
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

gene_freq_fun = function(y, freq) {
  # This function adds gene_freq and gene_length columns to mutations data frame
  freq$freq = with(freq, count / amino_acid_length)
  freq$freq = sapply(freq$freq, function(x) ifelse(is.nan(x),0,x))
  
  tmp_output = foreach (i=1:nrow(y), .combine = rbind) %dopar% {
    row_num = which(freq$gene_name == y$gene_name[i])
    if (length(row_num)) {
      c(freq$freq[row_num], freq$length[row_num])
    } else {
      c(NA, NA)
    }
  }
  y$gene_freq = tmp_output[,1]
  y$gene_length = tmp_output[,2]
  y = y[! is.na(y$gene_freq), ]
  
  return(y)
}

indel_filter = function(y, mega) {
  # flag indels close to meganormals
  l = dim(y)[1]
  y$mega = foreach (i = 1:l, .combine = c) %dopar% {
    mega_fun(y$chrom[i], y$pos[i], mega)
  }
  
  #   table(y$mega[y$indel == 1])
  #   table(y$mega[y$indel == 1 & y$Y == "som"])
  
  y$ref_len = nchar(y$ref)
  y$alt_len = nchar(y$alt)
  
  y$is_1 = as.numeric(y$is_1)
  y$is_2 = as.numeric(y$is_2)
  
  # only deletions
  y = y[ y$indel == "." | (y$ref_len > y$alt_len), ]
  
  # only Frame_shift
  y = y[ y$indel == "." | (y$effect == "FRAME_SHIFT"), ]
  
  # indels longer than 10 are excluded
  y = y[ y$ref_len <= 10 & y$alt_len <= 10, ]
  
  # remvoing indels close to meganormal
  y = y[ y$indel == "." | y$mega == 0, ]
  
  # removing indels with is_1 < 100
  y = y[ y$indel == "." | y$is_1 < 100, ]
  
  #   table(y$Y[y$indel == "1"])
  
  y$ref_len = NULL
  y$alt_len = NULL
  
  y
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

id2 = function(mt) {
  # .
  mt$in_none = as.numeric(! grepl("rs", mt$id) & ! grepl("COS", mt$id))*1   
  # only dbSNP
  mt$in_dbsnp = as.numeric(grepl("rs", mt$id) & ! grepl("COS", mt$id))*0
  # only Cosmic
  mt$in_cosmic = as.numeric(grepl("COS", mt$id) & ! grepl("rs", mt$id))*3
  # both Cosmic and dbSNP
  mt$in_both = as.numeric(grepl("rs", mt$id) & grepl("COS", mt$id))*2       
  
  mt$in_none + mt$in_dbsnp + mt$in_cosmic + mt$in_both
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

pfam_cosmic_fun = function(y, pfam_cosmic) {
  # This function adds the total number of cosmic mutations in a given domain
  # to the mutations that happen within that domain. The rest will have value of 0
  library(doMC)
  registerDoMC(cores = 4)
  
  y$pfam_cosmic = foreach (i = 1:dim(y)[1], .combine = "c") %dopar% {
    ind = (pfam_cosmic$chromStart[ pfam_cosmic$chrom == y$chrom[i]] <= y$pos[i]) &
      (pfam_cosmic$chromEnd[ pfam_cosmic$chrom == y$chrom[i]] >= y$pos[i])
    if (any(ind)) {
      max(pfam_cosmic$pval[ind])
    } else {
      0
    }
  }
  y
}

somatic_annot = function(y, somatic) {
  # somatic_file is a list of somatic mutations with "case_id", "chr", and
  # "start_position"
  
  y_list = unique(y$case)
  som_list = unique(somatic$case_id)
  
  # Choose the corresponding cases from somatic
  somatic = subset(somatic, somatic$case_id %in% y_list)
  
  somatic$ID = paste(somatic$case_id, somatic$chr, 
                     somatic$start_position, somatic$variant_allele, sep = ":")
  y$ID = paste(y$case, y$chrom, y$pos, y$alt, sep = ":")
  
  # Finding the somatic mutations using CBio data
  somatic$Not_Present = 1
  y$Y = "non_som"
  
  for (id in somatic$ID) {
    a = which(y$ID == id)
    if (length(a)) {
      y$Y[a] = "som"
      somatic$Not_Present[somatic$ID == id] = 0
    }
  }
  
  list(y, somatic)
}

tot_cosm_fun = function(y, cosm) {
  l = dim(y)[1]
  y$tot_cosm_samp = foreach (i = 1:l, .combine = c) %dopar% {
    row_num = which(cosm$gene == y$gene_name[i])
    if (length(row_num) > 0) {
      cosm$tot_nsample[row_num]
    } else {
      0
    }
  }
  y
}