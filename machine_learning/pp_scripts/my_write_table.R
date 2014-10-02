my_write_table = function(variable, output_file) {
    write.table(variable, file = output_file,
                quote = FALSE, 
                sep = "\t", 
                na = ".", 
                row.names = FALSE, 
                col.names = TRUE)
    cat("- File written.\n")
}