################################ ggStruct: structure visualization with ggplot2 ####################################
##run with < ggstruct("path to str output","path to str input",n.samples,vector of colors,k) >
##structure input file must include putative population origin in second column and have each individual on 2 rows
##plots will always be in the order of the structure data file (L to R for top to bottom)
##colors should be a vector of color names (strings) in the order of the desired colors for each population as 
##written in the putative pop origin column of the input file. 

ggstruct <- function(strOutput,strInput,n.samples,colors,k) {
  require(ggplot2);require(plyr);require(reshape);
  samples <- read.table(strInput,sep="\t")[1:2]
  samples <- samples[c(TRUE,FALSE),]
  colnames(samples) <- c("sample","population")
  s <- n.samples
  infile <- readLines(strOutput)
  a <- grep("Inferred ancestry of individuals:",infile)
  b <- grep("Estimated Ln Prob of Data",infile)
  lnml <- infile[b]
  df <- infile[(a+2):(a+s+1)]
  df <- strsplit(df,":  ")
  df <- do.call(rbind.data.frame, df)
  df <- df[2]
  colnames(df) <- c("ancestry")
  df$ancestry <- as.character(df$ancestry)
  df <- colsplit(df$ancestry,split=" ",names=c(1:k))
  df <- cbind(df,samples,row.names=NULL)
  df <- melt(df,id.var=c("sample","population"))
  popColors <- c()
  for (i in (1:k)){
    test <- subset(df,df$population==i)
    a <- test$variable[which(test$value == max(test$value))][1]
    popColors <- append(popColors,as.character(a))
  }
  palette<-c()
  palette[as.numeric(substr(popColors[1],2,3))] <- colors[1]
  palette[as.numeric(substr(popColors[2],2,3))] <- colors[2]
  palette[as.numeric(substr(popColors[3],2,3))] <- colors[3]
  palette[as.numeric(substr(popColors[4],2,3))] <- colors[4]
  palette[as.numeric(substr(popColors[5],2,3))] <- colors[5]
  palette[which(is.na(palette == TRUE))] <- colors[which(palette %in% colors == FALSE)]
  palette[duplicated(palette)] <- colors[which(colors %in% palette == FALSE)]
  
  ggplot(data=df,aes(x=sample,y=value,fill=variable,frame=population))+
    theme_bw()+theme(axis.text.x=element_text(angle=50,hjust=1.15,vjust=1.08,size=8),title=element_text(size=8),
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    scale_fill_manual(values=palette,name="Percent Ancestry")+xlab(NULL)+ylab(NULL)+
    geom_bar(stat="identity")+
    ggtitle(lnml)
  }

