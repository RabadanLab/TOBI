na_to_0 = function(x) {
  # This function converts a vector into numbers, and the resulting NAs into -1
  x = as.numeric(x)
  x[is.na(x)] = 0 
  x
}
