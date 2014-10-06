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