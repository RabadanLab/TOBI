id2 = function(mt) {
    # .
    mt$in_none = as.numeric(! grepl("rs", mt$id) & ! grepl("COS", mt$id))*1   
    # only dbSNP
    mt$in_dbsnp = as.numeric(grepl("rs", mt$id) & ! grepl("COS", mt$id))*0
    # only Cosmic
    mt$in_cosmic = as.numeric(grepl("COS", mt$id) & ! grepl("rs", mt$id))*3
    # both Cosmic and dbSNP
    mt$in_both = as.numeric(grepl("rs", mt$id) & grepl("COS", mt$id))*2       
    
    mt$in_none + mt$in_dbsnp + mt$in_cosmic + mt$in_both
}