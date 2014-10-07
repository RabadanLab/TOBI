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