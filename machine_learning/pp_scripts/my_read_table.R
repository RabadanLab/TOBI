my_read_table = function(x) {
    read.table(x, 
               sep = "\t", 
               header = TRUE, 
               stringsAsFactors = FALSE)
}