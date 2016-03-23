num_mut_gene = function(x) {
  # this function calculates the total number of mutation in a gene for all samples, 
  # and adds it to the mutation as a feature
  # x is the dataframe
  # output is a vector
  # suggested usage: x$num_mut_gene = num_mut_gene(x)
  gene_mut_freq = as.data.frame(table(x$gene_name))
  y = sapply(x$gene_name, 
             function(x) gene_mut_freq$Freq[gene_mut_freq$Var1 == x])
  num_cases = length(unique(x$case))
  y / num_cases * 100
}