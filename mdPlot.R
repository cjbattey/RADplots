###Plot % missing data in a phylip file###

mdPlot <- function(file,arrange=T){
    require(data.table);require(stringr);require(plyr)
  a <- fread(file,header=T)
  nchar <- as.numeric(names(a)[2])
  colnames(a) <- c("sample","seq")
  b <- lapply(a$seq,FUN=function(e) str_count(e,"N"))
  c <- as.numeric(b)/nchar
  d <- data.frame(sample=a$sample,md=c)
  if(arrange==F){
    barplot(d$md,names.arg=d$sample,las=2,cex.names=0.75,ylim=c(0,1))
  }
    d <- arrange(d,md)
    barplot(d$md,names.arg=d$sample,las=2,cex.names=0.75,ylim=c(0,1))
}
