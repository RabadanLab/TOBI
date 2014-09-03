rm(list = ls())

################################################################################################################################################
### functions
################################################################################################################################################
is.letter <- function(x) grepl("[[:alpha:]]", x)

##################################################################################################################
comma_to_num = function(x) {
  y = sapply(x, function(x) {
    max(as.numeric(unlist(strsplit(as.character(x), ","))))
  })
  return(unname(y))
}

################################################################################################################################################
### pre-processing
### adding the Y label
################################################################################################################################################

normal = 0

somatic = read.table("/Volumes/ifs/scratch/Results/Somatic_Mutations/CBio/cbio_somatic_SNPs.txt", 
                     sep = "\t", 
                     header = TRUE, 
                     stringsAsFactors = FALSE)
gene_freq_length = read.table("/Volumes/ifs/scratch/meganormal/gene_count_acl_length_merged.txt", 
                              sep = "\t", 
                              header = TRUE, 
                              stringsAsFactors = FALSE)
cbio_cases = readLines("/Volumes/ifs/scratch/Results/Somatic_Mutations/CBio/All_cases.txt")

mt = read.table("/Volumes/ifs/scratch/Results/GBM/final_tables/all_GBM_mutations_Aug-6-2014_filt_indel_techn_biol_no_meg2.txt", 
                sep = "\t", 
                header = TRUE, 
                stringsAsFactors = FALSE)


cnv = read.table("/Volumes/ifs/scratch/Results/Somatic_Mutations/CBio/cbio_cna.txt", 
                 sep = "\t", 
                 header = TRUE, 
                 stringsAsFactors = FALSE)
cnv$GENE_ID = NULL
colnames(cnv) = gsub("\\.", "-", colnames(cnv))

tot_cosm_samp = read.table("/Volumes/ifs/scratch/cosmic/tot_num_samples_gene.txt", 
                           sep = "\t", 
                           header = TRUE, 
                           stringsAsFactors = FALSE)



# Exception list
exc_list = NULL
# exc_list = c("TCGA-06-0132", "TCGA-06-0178", "TCGA-06-0165", "TCGA-19-5952", 
#              "TCGA-06-0167", "TCGA-06-5410", "TCGA-32-1980", "TCGA-06-0141", "TCGA-06-0189")
mt_list = unique(mt$case)
mt_list = mt_list[! mt_list %in% exc_list]
mt_list = mt_list[mt_list %in% cbio_cases]

# Choose the corresponding mt from somatic
somatic = subset(somatic, somatic$case_id %in% mt_list)
mt = subset(mt, mt$case %in% mt_list)

somatic$ID = paste(somatic$case_id, somatic$chr, somatic$start_position, sep = ":")
mt$ID = paste(mt$case, mt$chrom, mt$pos, sep = ":")

# Finding the somatic mutations using CBio data
somatic$Not_Present = 0
mt$Y = 0

for (id in somatic$ID) {
  if (length(which(mt$ID == id))) {
    mt$Y[which(mt$ID == id)] = 1
  } else {
    somatic$Not_Present[which(somatic$ID == id)] = 1
  }
}

# Removing redundant columns from cnv
cnv = cnv[, colnames(cnv) %in% c("COMMON", mt_list)]

# checking which mutations in important genes are missed by filtering
list_driver_genes = c("PIK3R1", "PTEN", "PIK3CA", "TP53", "EGFR", "IDH1", "BRAF", "RB1", "NF1", "PDGFRA", "LRP2", "PPP2R3A", "FAT2",
                      "GPR116", "RIMS2", "NOS1", "PIK3C2B", "MDM4", "MYCN", "KIT", "DDIT3", "TSHZ2", "LRP1B", "CTNND2", "HCN1", "PKHD1",
                      "TEK", "PCNX", "HERC2", "LZTR1", "BCOR", "ATRX", "PCDH11X", "CDK4", "CDKN2A")

missed_somatic_genes = unique(somatic$gene_symbol[somatic$Not_Present == 1])

(missed_driver_genes = missed_somatic_genes[missed_somatic_genes %in% list_driver_genes])

View(somatic[somatic$Not_Present == 1 & somatic$gene_symbol %in% missed_driver_genes, ])

table(mt$Y)
table(somatic$Not_Present)


# writing the mutations to be checked
missed_driver = somatic[somatic$Not_Present == 1 & somatic$gene_symbol %in% missed_driver_genes, ]
missed_driver = missed_driver[, c("gene_symbol", "case_id", "chr", "start_position")]
write.table(missed_driver, file = "/Volumes/ifs/scratch/Results/Finding_lost_mutations/missed_driver_mutations_lib.txt", quote = FALSE, sep = "\t", na = ".", row.names = FALSE, col.names = FALSE)


# adding copy number variantion, and total number of cosmic sample data
l = dim(mt)[1]
mt$cnv = NA
mt$tot_cosm_samp = 0
if (normal == 1) {
  mt$cnv = 1
} else {
  for (i in 1:l) {
    # adding cnv
    row_num = which(mt$gene_name[i] == cnv$COMMON)
    col_num = which(mt$case[i] == colnames(cnv))
    if (length(row_num) == 1) {
      if (is.nan(cnv[row_num, col_num])) {
        mt$cnv[i] = NA
      } else{
        mt$cnv[i] = 2 ^ cnv[row_num, col_num]
      }
    } else {
      mt$cnv[i] = NA
    }
  }
}

for (i in 1:l) {
  # adding total number of cosmic samples to the data
  row_num = which(tot_cosm_samp$gene == mt$gene_name[i])
  if (length(row_num) > 0) {
    mt$tot_cosm_samp[i] = tot_cosm_samp$tot_nsample[row_num]
  }
}

# Calculating mutation allele frequency
mt$freq = (mt$dp4_3 + mt$dp4_4) / (mt$dp4_1 + mt$dp4_2 + mt$dp4_3 + mt$dp4_4) * 100

# Correcting factors

mt$effect_impact[mt$effect_impact == "."] = 0
mt$effect_impact[mt$effect_impact == "LOW"] = 1
mt$effect_impact[mt$effect_impact == "MODERATE"] = 2
mt$effect_impact[mt$effect_impact == "MODIFIER"] = 3
mt$effect_impact[mt$effect_impact == "HIGH"] = 4

mt$gt_1[mt$gt_1 == '0/1'] = 0
mt$gt_1[mt$gt_1 == '1/1'] = 1

# Removing empty columns

mt$indel = NULL
mt$is_1 = NULL
mt$is_2 = NULL
mt$common = NULL
mt$g5 = NULL
mt$g5a = NULL
mt$nov = NULL
mt$meganormal_id = NULL
# mt$src = NULL
mt$vp = NULL
mt$transcript_id = NULL
mt$mut = NULL
mt$mtp = NULL
mt$nsf = NULL
mt$tpa = NULL
mt$kgvalidated = NULL
mt$oth = NULL
mt$nmutperid = NULL
mt$dbNSFP_1000Gp1_AC = NULL
mt$dbNSFP_1000Gp1_AFR_AC = NULL
mt$dbNSFP_1000Gp1_AMR_AC = NULL
mt$dbNSFP_1000Gp1_ASN_AC = NULL
mt$dbNSFP_1000Gp1_EUR_AC = NULL
mt$dbNSFP_CADD_raw = NULL
mt$dbNSFP_CADD_raw_rankscore = NULL
mt$dbNSFP_UniSNP_ids = NULL
mt$dbSNPBuildID = NULL
mt$filter = NULL
mt$an = NULL
mt$aa = NULL
mt$ac = NULL
mt$af1 = NULL
mt$asp = NULL
mt$ass = NULL
mt$caf = NULL
mt$cda = NULL
mt$cds = NULL
mt$cfl = NULL
mt$cgt = NULL
mt$clnacc = NULL
mt$clndsdb = NULL
mt$clndsdbid = NULL
mt$clnhgvs = NULL
mt$clnsrc = NULL
mt$clnsrcid = NULL
mt$clr = NULL
mt$dss = NULL
mt$genotype_number = NULL
mt$g3 = NULL
mt$geneinfo = NULL
mt$gene = NULL
mt$hwe = NULL
mt$indel.1 = NULL
mt$lof = NULL
mt$mdv = NULL
mt$nmd = NULL
mt$noc = NULL
mt$nsn = NULL
mt$om = NULL
mt$pc2 = NULL
mt$pchi2 = NULL
mt$pmc = NULL
mt$pr = NULL
mt$qbd = NULL
mt$qchi2 = NULL
mt$r3 = NULL
mt$rs = NULL
mt$rspos = NULL
mt$strand = NULL
mt$u3 = NULL
mt$u5 = NULL
mt$ugt = NULL
mt$wtd = NULL
mt = mt[, -grep("*rankscore*", colnames(mt))]
mt = mt[, -grep("*pred*", colnames(mt))]
mt = mt[, -grep("*cln*", colnames(mt))]

mt$sao[mt$sao == "."] = "NA"
mt$sao[mt$sao == "0"] = "Unspecified"
mt$sao[mt$sao == "1"] = "Germline"

mt$warnings[mt$warnings != "."] = 1
mt$warnings[mt$warnings == "."] = 0

mt$som_counts = as.numeric(as.character(mt$som_counts))
mt$rpb = as.numeric(as.character(mt$rpb))
mt$exon_rank = as.numeric(as.character(mt$exon_rank))
mt$amino_acid_length = as.numeric(as.character(mt$amino_acid_length))
# mt$cosmic_nsamp = as.numeric(as.character(mt$cosmic_nsamp))
mt$ssr = as.numeric(as.character(mt$ssr))
mt$fq = as.numeric(mt$fq)

gene_freq_length$freq = gene_freq_length$count / gene_freq_length$amino_acid_length
gene_freq_length$freq = sapply(gene_freq_length$freq, function(x) ifelse(is.nan(x),0,x))

mt$gene_freq = 0
mt$gene_length = 0

for (i in 1:dim(mt)[1]) {
  row_num = which(gene_freq_length$gene_name == mt$gene_name[i])
  if (length(row_num)) {
    mt$gene_freq[i] = gene_freq_length$freq[row_num]
    mt$gene_length[i] = gene_freq_length$length[row_num]
  }
}

# removing the abnormals
if (normal != 1) {
  num_cases = length(table(mt$case))
  y = data.frame(case = names(table(mt$case)), non_som = rep(0, num_cases), som = rep(0, num_cases))
  
  for (case in unique(mt$case)) {
    y$non_som[y$case == case] = table(mt$Y[mt$case == case])[1]
    y$som[y$case == case] = table(mt$Y[mt$case == case])[2]
  }
  
  exclude_list = as.character(y$case[is.na(y$som)])
  flag = ! mt$case %in% exclude_list
#   mt = mt[flag,]
}


mt$dbNSFP_CADD_phred = comma_to_num(mt$dbNSFP_CADD_phred)
mt$dbNSFP_1000Gp1_ASN_AF = NULL
mt$dbNSFP_ESP6500_AA_AF = comma_to_num(mt$dbNSFP_ESP6500_AA_AF)
mt$dbNSFP_phastCons100way_vertebrate = comma_to_num(mt$dbNSFP_phastCons100way_vertebrate)
mt$dbNSFP_1000Gp1_EUR_AF = comma_to_num(mt$dbNSFP_1000Gp1_EUR_AF)
mt$dbNSFP_1000Gp1_AFR_AF = NULL
mt$dbNSFP_ESP6500_EA_AF = comma_to_num(mt$dbNSFP_ESP6500_EA_AF)
mt$dbNSFP_1000Gp1_AMR_AF = NULL
mt$dbNSFP_1000Gp1_AF = NULL
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

mt$ID = NULL
mt$pl_1 = NULL
mt$dbNSFP_SiPhy_29way_pi = NULL

# calculating the pvalue for being somatic based on dp and freq
for (i in 1: dim(mt)[1]) {
  var_dp = round(mt$dp[i]*mt$freq[i]/100)
  mt$pval_som[i] = dbinom(var_dp, mt$dp[i], 0.5)
}

# mt$logit = log10(mt$pval_som/(1-mt$pval_som))

mt$dbNSFP_1000Gp1 = as.numeric(! is.na(mt$dbNSFP_1000Gp1_EUR_AF))
mt$dbNSFP_1000Gp1_EUR_AF = NULL

# making id2
mt$in_none = as.numeric(! grepl("rs", mt$id) & ! grepl("COS", mt$id))*1   # .
mt$in_dbsnp = as.numeric(grepl("rs", mt$id) & ! grepl("COS", mt$id))*0    # only dbSNP
mt$in_cosmic = as.numeric(grepl("COS", mt$id) & ! grepl("rs", mt$id))*3   # only Cosmic
mt$in_both = as.numeric(grepl("rs", mt$id) & grepl("COS", mt$id))*2       # both Cosmic and dbSNP
mt$id2 = mt$in_none + mt$in_dbsnp + mt$in_cosmic + mt$in_both
mt$in_none = NULL
mt$in_dbsnp = NULL
mt$in_cosmic = NULL
mt$in_both = NULL

mt$Y[mt$Y == "0"] = "non_som"
mt$Y[mt$Y == "1"] = "som"

write.table(mt, file = "/Volumes/ifs/scratch/Results/GBM/final_tables/all_GBM_mutations_104cases_filt_indel_techn_biol_pre_proc.txt", 
            quote = FALSE, sep = "\t", na = ".", row.names = FALSE)

