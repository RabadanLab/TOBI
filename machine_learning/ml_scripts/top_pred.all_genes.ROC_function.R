driver_mutations_gene = function(mydata, somatic_path, out_path, output_plot = FALSE) {
  # this function counts genes that are predicted somatic, and gives the numbers in true somatic set, then TP, FP, FN,
  # and plots the outcome per gene
  # mydata is a dataframe
  library(dplyr)
  require(ggplot2)
  require(reshape2)

#  args <- commandArgs(TRUE)
#  variant_path <- args[[1]] # all variants after TOBI filters and predictions
#  somatic_path <-  args[[2]] # path to true somatic calls
#  out_path <- args[[3]] # path to final file
#  print(args)
#  
#  mydata <- read.table(variant_path, header = T, stringsAsFactors = F)
 
  # somatic mutations
  somatic = read.table(somatic_path, 
                       sep = "\t", 
                       header = TRUE, 
                       stringsAsFactors = FALSE)
  
  # Choose the corresponding tumor from somatic
  mydata_list = unique(mydata$case)
  somatic = subset(somatic, somatic$case_id %in% mydata_list)

  # generate performance values
  time1 = Sys.time()
  #from true somatic calls
  som_notes <- somatic %>% group_by(gene_symbol) %>% mutate(tot_true_som= n()) %>% 
    rename(gene_name=gene_symbol) %>% select(gene_name, tot_true_som) %>% distinct()
  #true calls after TOBI filtering
  true_after_filt <- mydata %>% filter(Y== "som") %>% group_by(gene_name, Y) %>% mutate(true_after_filt= n()) %>%  
    ungroup() %>%  select(gene_name, true_after_filt) %>% distinct()
  #predicted as somatic from TOBI
  tot_pred_som <- mydata %>% filter(pred== "som") %>% group_by(gene_name, pred) %>% mutate(tot_pred_som= n()) %>%  
    ungroup() %>%  select(gene_name, tot_pred_som) %>% distinct()
  #number of unique cases with true somatic mutations in the gene
  true_case <- mydata %>% filter(Y== "som") %>% select(gene_name, case, Y) %>% distinct() %>% 
    group_by(gene_name) %>% mutate(true_case= n()) %>%  
    ungroup() %>%  select(gene_name, true_case) %>% distinct()
  #number of unique cases with predicted somatic mutations in the gene
  pred_case <- mydata %>% filter(pred== "som") %>% select(gene_name, case, pred) %>% distinct() %>% 
    group_by(gene_name) %>% mutate(pred_case= n()) %>%  
    ungroup() %>%  select(gene_name, pred_case) %>% distinct()
  TP <- mydata %>% filter(Y== "som" & pred== "som") %>% group_by(gene_name) %>% mutate(TP= n()) %>% 
    select(gene_name, TP) %>% distinct()
  FP <- mydata %>% filter(Y== "non_som" & pred== "som") %>% group_by(gene_name) %>% mutate(FP= n()) %>% 
    select(gene_name, FP) %>% distinct()
  FN <- mydata %>% filter(Y== "som" & pred== "non_som") %>% group_by(gene_name) %>% mutate(FN= n()) %>% 
    select(gene_name, FN) %>% distinct()
  
  #select longest transcript for protein length since I've lost transcript labelling data
  protein_length <- mydata %>% rename(protein_length=amino_acid_length) %>% 
    select(gene_name, protein_length) %>% distinct() %>% group_by(gene_name) %>% top_n(1,)  

    # fraction of synonymous mutations over all
  syn<- mydata %>% filter(effect == "SYNONYMOUS_CODING" | effect == "SYNONYMOUS_STOP") %>% group_by(gene_name) %>% 
    mutate(syn=n()) %>% select(gene_name,syn) %>% distinct()
  tot_length <- mydata %>% group_by(gene_name) %>% mutate(tot_length=n()) %>% select(gene_name,tot_length) %>% distinct()
  syn_ratio <- full_join(syn, tot_length, "gene_name") %>% mutate_each(funs(replace(., which(is.na(.)), 0))) %>% 
    mutate(syn_ratio=syn/tot_length) %>% select(gene_name,syn_ratio) %>% distinct()

# all cases with mutations after filtering
case_gene <- mydata %>% select(gene_name,case) %>% distinct() %>% group_by(gene_name) %>% 
  mutate(all_cases=n()) %>% select(gene_name,all_cases) %>% distinct()

# join all above categories 
  imp_mt <-  right_join(full_join(full_join(full_join(full_join(full_join(full_join(full_join(full_join(full_join(som_notes, true_after_filt, "gene_name"),
                      tot_pred_som, "gene_name"), 
                      true_case, "gene_name"), 
                      pred_case, "gene_name"),
                      TP, "gene_name"),
                      FP, "gene_name"), 
                      FN, "gene_name"), 
                      protein_length, "gene_name"),
                      syn_ratio, "gene_name"), 
                      case_gene, "gene_name")  %>%  
    mutate_each(funs(replace(., which(is.na(.)), 0)))
  print(  Sys.time() - time1)
  
  cat("Total number of predicted somatically altered variants:", sum(imp_mt$tot_pred_som), "\n")
  cat("Total true somatically altered variants:", sum(imp_mt$tot_true_som), "\n")
  cat("After filtering and variants:", sum(imp_mt$TP), "\n")

  write.table(imp_mt, out_path, 
              sep = "\t", 
              quote = FALSE, 
              na = ".", 
              row.names = FALSE)
  
  if (output_plot) {
    imp_mt_melt = melt(imp_mt, id = c("gene"))
    q = qplot(x = gene, y = value , fill=variable,
              data = imp_mt_melt, geom="bar", stat="identity",
              position="dodge")
    #   q
    q = q + 
      theme_bw() + 
      coord_fixed(ratio = 0.6) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), 
            text = element_text(size = 18),
            legend.title = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            legend.position="bottom")
    # q
    print(q)
  }
}