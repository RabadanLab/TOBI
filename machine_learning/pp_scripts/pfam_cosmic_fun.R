pfam_cosmic_fun = function(y, pfam_cosmic) {
  # This function adds the total number of cosmic mutations in a given domain
  # to the mutations that happen within that domain. The rest will have value of 0
  library(doMC)
  registerDoMC(cores = 4)
  
  y$pfam_cosmic = foreach (i = 1:dim(y)[1], .combine = "c") %dopar% {
    ind = (pfam_cosmic$chromStart[ pfam_cosmic$chrom == y$chrom[i]] <= y$pos[i]) &
      (pfam_cosmic$chromEnd[ pfam_cosmic$chrom == y$chrom[i]] >= y$pos[i])
    if (any(ind)) {
      max(pfam_cosmic$all_cosmic[ind])
    } else {
      0
    }
  }
  y
}