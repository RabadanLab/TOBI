comma_to_num = function(x) {
  # This function cleans a vector with multiple values in random cells,
  # and calls na_minus1 function to replace the NAs with -1
  y = sapply(x, function(x) {
    max(as.numeric(unlist(strsplit(as.character(x), ","))))
  })
  y = unname(y)
  y = na_minus1(y)
  return(unname(y))
}