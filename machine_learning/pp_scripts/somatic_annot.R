somatic_annot = function(y, somatic) {
  # somatic_file is a list of somatic mutations with "case_id", "chr", and
  # "start_position"
  
  y_list = unique(y$case)
  som_list = unique(somatic$case_id)
  
  # Choose the corresponding cases from somatic
  somatic = subset(somatic, somatic$case_id %in% y_list)
  
  somatic$ID = paste(somatic$case_id, somatic$chr, somatic$start_position, sep = ":")
  y$ID = paste(y$case, y$chrom, y$pos, sep = ":")
  
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