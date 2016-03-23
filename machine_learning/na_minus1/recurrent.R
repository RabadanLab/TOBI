recurrent = function(y) {
  require(dplyr)
  
  y %>% 
    mutate(uniq = paste(chrom, pos, ref, alt, sep = ":")) -> y
  
  y %>%
    group_by(uniq) %>%
    summarize(recur = n()) %>%
    inner_join(y, ., by = "uniq") %>%
    select(- uniq) -> y
  
  return(y$recur)
}