getproteinlabeldict <- function(){ 
con <- file("./www/data/GESA_Canonical_pathways_c2.cp.v5.0.symbols.gmt_filtered.tsv") 
open(con);
results.list <- list();
current.line <- 1
map <- new.env(hash=T, parent=emptyenv())
while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
  
  dt <- strsplit(line, split="\t")
  val<-unlist(dt[[1]][1])
  keys<-unlist(dt[[1]][-1])
  for(key in keys){
    map[[key]] <- val
  }
  #results.list[[current.line]] <- toString(unlist(strsplit(line, split="\t")))
  current.line <- current.line + 1
} 
close(con)

return(map)
}