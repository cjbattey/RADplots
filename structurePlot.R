#################################################################################################
##################### structurePlot: lightweight structure output plotter #######################
#################################################################################################
##Required input: structure output file
##options: colors = vector of at least k colors. Default has max 8 colors. 
##         pop.info = TRUE/FALSE for presence of column of putative pop info. Defaults to TRUE, required if FALSE. 
##         pop.order = FALSE (default) plots in order of structure input file
##         pop.order = TRUE plots in order of populations as in structure input file optional column "putative pop info"
#Note: brewer.pal will throw a warning on k=2 runs. Don't worry about it - it's just telling you there's an extra color in the default palette.
#Please cite github URL if you use this in a publication

structurePlot <- function(strOutput,pop.info=TRUE,pop.order=FALSE,colors=brewer.pal(k,name="Set1")) {
  require(plyr);require(reshape);require(RColorBrewer)
  infile <- readLines(strOutput,warn=FALSE)
  a <- grep("Inferred ancestry of individuals:",infile)+2
  b <- grep("Estimated Ln Prob of Data",infile)
  c <- grep("Estimated Allele Frequencies in each cluster",infile)-3
  k <- as.numeric(substr(infile[grep("populations assumed",infile)],4,5))
  lnml <- infile[b]
  df <- infile[(a):(c)]
  df <- strsplit(df," +")
  df <- do.call(rbind.data.frame, df)
  if (pop.info=="TRUE" & pop.order=="FALSE"){
    df <- data.frame(df[3],df[5],df[7:(6+k)])
    colnames(df) <- c("sample","pop",1:k)
    barplot(t(as.matrix(df[3:(2+k)])),col=colors,border=NA,names=df$sample,cex.names=0.75,las=2,main=lnml,cex.main=0.75,font.main=1)
  }
  if (pop.info=="TRUE" & pop.order=="TRUE"){
    df <- data.frame(df[3],df[5],df[7:(6+k)])
    colnames(df) <- c("sample","pop",1:k)
    df <- arrange(df,df$pop)
    barplot(t(as.matrix(df[3:(2+k)])),col=colors,border=NA,names=df$sample,cex.names=0.75,las=2,main=lnml,cex.main=0.75,font.main=1)
  }
  if (pop.info=="FALSE"){
    df <- data.frame(df[3],df[6:(5+k)])
    colnames(df) <- c("sample",1:k)
    names <- df$sample
    sample <- "sample"
    df <- df[,!(names(df)) %in% sample]
    barplot(t(as.matrix(df)),col=colors,border=NA,names=names,cex.names=0.75,las=2,main=lnml,cex.main=0.75,font.main=1)
  } 
}


