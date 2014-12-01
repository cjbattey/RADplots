### plot missing data by sample across an alignment, as well as total missing data ###
### Input is a phylip formatted alignment like pyRAD's .phy and .unlinked_snps output, but NOT the .snp output b/c of the extra spaces bw each loci's SNP's
### Note: this function is pretty slow (~2 min for 20,000 SNP alignment). Would be way faster without the for loop. 

mdPlot <- function(file) {
  require(ggplot2);require(reshape);require(plyr)
  data <- read.table(file,header=TRUE)
  colnames(data) <- c("sample","sequence")
  data$sequence <- as.character(data$sequence)
  data <- data.frame(colsplit(data$sequence,split="",names=c(1:nchar(data$sequence[1]))),row.names=data$sample)
  samples <- c(rownames(data))
  missingData <- c()
  for(i in 1:nrow(data)) {
    a <- data[i,]
    a <- a[which(a == "N" | a == "?" | a == "-")]
    b <- length(a)
    missingData <- append(missingData,b)
  }
  md <- data.frame(samples,missingData)
  md$missingData <- 100*(md$missingData/length(data))
  md <- arrange(md,desc(md$missingData))
  md$samples <- reorder(md$samples,desc(md$missingData))
  totalMD <- 100*(sum(missingData)/(length(data)*nrow(data)))
  ggplot(data=md,aes(x=samples,y=missingData))+theme_bw()+xlab(NULL)+ylab("% Missing Data")+
    theme(axis.text.x=element_text(angle=90,hjust=1.0,vjust=0.5))+
    geom_bar(stat="identity")+geom_hline(y=totalMD,col="red")+
    annotate(geom="text",label="Total Missing Data",col="red",y=totalMD+0.08,x=0.98*length(md$samples),cex=4)+coord_flip()
}
