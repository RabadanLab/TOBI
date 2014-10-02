case_mut_gene = function(x) {
  # this function calculates the total number of mutation in a gene for every case, 
  # and adds it to the mutation as a feature
  # x is the dataframe
  # output is a vector
  # suggested usage: x$case_mut_gene = case_mut_gene(x)
  x$case_mut_gene = NA
  for (item in unique(x$case)) {
    gene_mut_freq = as.data.frame(table(x$gene_name[x$case == item]))
    x$case_mut_gene[x$case == item] = sapply(x$gene_name[x$case == item], 
                                             function(x) gene_mut_freq$Freq[gene_mut_freq$Var1 == x])
  }
  num_cases = length(unique(x$case))
  x$case_mut_gene / num_cases * 100
}