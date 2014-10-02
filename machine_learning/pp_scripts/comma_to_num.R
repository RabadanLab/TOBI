comma_to_num = function(x) {
    y = sapply(x, function(x) {
        max(as.numeric(unlist(strsplit(as.character(x), ","))))
    })
    return(unname(y))
}