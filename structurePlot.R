########## structurePlot: an R function for plotting Structure output
##Required input: structure output file
##options: colors = vector of at least k colors. Default has max 8 colors. 
##         pop.order = FALSE (default) plots in order of structure input file
##         pop.order = TRUE plots in order of populations as in structure input file optional column "putative pop info"
##         pop.info = TRUE/FALSE for presence of column of putative pop info. Defaults to TRUE, required if FALSE. 
##         spacing = Amount of white space left bw bars, as a fraction of bar width. Defaults to 0. R doesn't seem to do well with less than 0.2. 
##         outline = color used for outline of every bar (default is black). For no outline set to NA. 
##
##experimental: color.matching = TRUE will attach colors to putative populations (requires pop.info=TRUE) for consistent 
##plotting across multiple runs. The cluster attached to a given putative population is taken to be that with the highest 
##average ancestry across all individuals in that population. Works with my data but there are probably still bugs in here. 
##If color.matching=TRUE, user-specified colors are assigned by population # - first color in the vector for pop 1, second for pop 2, etc.
##
##For additional control over plotting (line width, spacing, etc), see ?barplot and modify the relevant barplot call in the function

structurePlot <- function(strOutput,pop.order=FALSE,pop.info=TRUE,colors=c("gold","forestgreen","magenta3","orangered","cornflowerblue","orange","sienna","dodgerblue4"),color.matching=FALSE,spacing=0,outline='black') {
  require(plyr);require(reshape)
  infile <- readLines(strOutput,warn=FALSE)
  a <- grep("Inferred ancestry of individuals:",infile)+2
  b <- grep("Estimated Ln Prob of Data",infile)
  c <- grep("Estimated Allele Frequencies in each cluster",infile)-3
  k <- as.numeric(substr(infile[grep("populations assumed",infile)],4,5))
  lnml <- infile[b]
  df <- infile[(a):(c)]
  df <- strsplit(df," +")
  df <- do.call(rbind.data.frame, df)
  if (pop.info=="TRUE"){
    df <- data.frame(df[3],df[5],df[7:(6+k)])
    colnames(df) <- c("sample","pop",1:k)
    n.pops <- max(as.numeric(as.character(df$pop)))
  if (color.matching=="TRUE"){
    palette <- c(rep("NA",k))
    for(i in 1:n.pops){
      d <- subset(df,pop==i)[3:(2+k)]
      for(j in 1:ncol(d)) {
        d[,j] <- as.numeric(as.character(d[,j]))
      }
      e <- colMeans(d)
      f <- as.numeric(names(e[which(e==max(e))]))
      if (palette[f] == "NA"){
        palette[f] <- colors[i]
        }
      }
    if (n.pops < k) {
      for (i in 1:length(colors)){
          palette[which(palette == "NA")] <- colors[which(colors%in%palette == FALSE)]
        }
      }
    colors <- palette
   }
  if (pop.order=="TRUE"){
    df <- arrange(df,df$pop)
    barplot(t(as.matrix(df[3:(2+k)])),axes=FALSE,col=colors,border=outline,names=df$sample,cex.names=0.75,las=2,main=lnml,cex.main=0.75,font.main=1,space=spacing,xpd=FALSE)
  }
    barplot(t(as.matrix(df[3:(2+k)])),axes=FALSE,col=colors,border=outline,names=df$sample,cex.names=0.75,las=2,main=lnml,cex.main=0.75,font.main=1,space=spacing,xpd=FALSE)
  if (pop.info=="FALSE"){
    df <- data.frame(df[3],df[6:(5+k)])
    colnames(df) <- c("sample",1:k)
    n.pops <- k
    names <- df$sample
    sample <- "sample"
    df <- df[,!(names(df)) %in% sample]
    barplot(t(as.matrix(df)),axes=FALSE,col=colors,border=outline,names=names,cex.names=0.75,las=2,main=lnml,cex.main=0.75,font.main=1,space=spacing,xpd=FALSE)
  }
  } 
}


