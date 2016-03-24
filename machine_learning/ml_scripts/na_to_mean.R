na_to_mean = function(x) {
  # This function converts a vector into numbers, 
  # and the resulting NAs into the mean of the present values 
  x = as.numeric(x)
  y = as.numeric(x[x!=-1])
  new = mean(y, na.rm = TRUE)
  x[x==-1] = new #-1
  x
}
