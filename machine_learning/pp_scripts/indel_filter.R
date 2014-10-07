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