driver_mutations = function(mydata, modboost, somatic_path, output_plot = TRUE) {
  # this function finds how many driver mutation are left after filtering and ML,
  # and plots the outcome per gene
  # mydata is a dataframe
  
  require(ggplot2)
  require(reshape2)
  
  mydata$pred = unlist(predict(modboost[[1]], mydata[,-c(1:11)]))
#   mydata$pred = predict(modboost, mydata[,-c(1:11)])
  
  GBM_Recur_Mutations = c("PIK3R1", "PTEN", "PIK3CA", "TP53", "EGFR", "IDH1", 
                          "BRAF", "RB1", "NF1", "PDGFRA", "LRP2", "PPP2R3A", 
                          "FAT2", "GPR116", "RIMS2", "NOS1", "PIK3C2B", "MDM4", 
                          "MYCN", "KIT", "DDIT3", "TSHZ2", "LRP1B", "CTNND2", 
                          "HCN1", "PKHD1", "TEK", "PCNX", "HERC2", "LZTR1", 
                          "BCOR", "ATRX", "PCDH11X", "CDK4", "CDKN2A")
  
  GBM_Recur_Mutations = GBM_Recur_Mutations[GBM_Recur_Mutations %in% 
                                              unique(mydata$gene_name)]
  
  # somatic mutations
  somatic = read.table(somatic_path, 
                       sep = "\t", 
                       header = TRUE, 
                       stringsAsFactors = FALSE)
  
  # Choose the corresponding tumor from somatic
  mydata_list = unique(mydata$case)
  somatic = subset(somatic, somatic$case_id %in% mydata_list)
  
  l = length(GBM_Recur_Mutations)
  imp_mt = data.frame(gene = GBM_Recur_Mutations, 
                      tot_num_mut = rep(0, l), 
                      num_mut_filt = rep(0, l), 
                      TP = rep(0, l), 
                      FP = rep(0, l), 
                      FN = rep(0, l))
  for (i in 1:l) {
    imp_mt$tot_num_mut[i] = sum(somatic$gene_symbol == imp_mt$gene[i])
    imp_mt$num_mut_filt[i] = sum(mydata$gene_name == imp_mt$gene[i] & 
                                   mydata$Y == "som")
    imp_mt$TP[i] = sum(mydata$gene_name == imp_mt$gene[i] & 
                         mydata$Y == "som" & 
                         mydata$pred == "som")
    imp_mt$FP[i] = sum(mydata$gene_name == imp_mt$gene[i] & 
                         mydata$Y == "non_som" & 
                         mydata$pred == "som")
    imp_mt$FN[i] = sum(mydata$gene_name == imp_mt$gene[i] & 
                         mydata$Y == "som" & 
                         mydata$pred == "non_som")
  }
  imp_mt = imp_mt[imp_mt$tot_num_mut != 0, ]
  cat("Total number of driver mutations:", sum(imp_mt$tot_num_mut), "\n")
  cat("After filter:", sum(imp_mt$num_mut_filt), "\n")
  cat("After ML:", sum(imp_mt$TP), "\n")
  
  if (output_plot) {
    imp_mt_melt = melt(imp_mt, id = c("gene"))
    q = qplot(x = gene, y = value, fill=variable,
              data = imp_mt_melt, geom="bar", stat="identity",
              position="dodge")
    q = q + 
      theme_bw() + 
      coord_fixed(ratio = 0.6) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), 
                               text = element_text(size = 18),
                               legend.title = element_blank(),
                               axis.title.x = element_blank(),
                               axis.title.y = element_blank(),
                               legend.position="bottom")
    print(q)
  }
}