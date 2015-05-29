### Write D test input files in pyRAD format ###
#list permutations of ingroup taxa, with gilvus as the outgroup
#skip if p1=p2
#skip if [p2] [p1] [p3] is already present (IE no swapping ingroups)
#need to add extra formatting commands block at top - see old files or run ../pyrad -D
library(reshape)

olivN<-c("olivN_BTS07238_MX","olivN_JK01431_FL","olivN_AWH351_AZ","olivN_JK03328_CAN","olivN_JK03139_NY","olivN_TNL197_TX","olivN_DHB5712_LOU","olivN_JMD970_CA")
olivS<-c("olivS_69665_TTO","olivS_69311_TTO","olivS_69419_TTO","olivS_627389_GUY","olivS_392424_BRZ","olivS_GAV834_ARG","olivS_334603_BOL","olivS_427316_BRA","olivS_427315_BRA2","olivS_632655_GUY")
alt<-c("alt_13307_LOU","alt_21786_LOU","alt_331141_JAM","alt_3814_LOU","alt_DOT6379_DR","alt_DOT6385_DR")
mag<-c("mag_DOT18067_MX","mag_DOT18078_MX","mag_DOT18080_MX")
flav<-c("flav_BTS08151_MXW","flav_BTS08203_MXW","flav_BTS08344_MXE","flav_BTS08346_MXE","flav_DHB3401_MXE","flav_EAG037_MXW","flav_JK02216_GUA","flav_SAR7936_MXW")
gilv<-c("gilv_GMS898_MX","gilv_JK03429_CA")
samples <- list(olivN,olivS,alt,mag,flav,gilv)

command <- c()
for(b in 1:5){
  for(c in 1:5)
    if (c != b){
      for (i in samples[[b]]){
        for(j in samples[[b]]){
          for(k in samples[[c]]){
            #skip final loop if using fixed combo outgroup
              p <- grep(paste("\\[",j,"\\]","\t","\\[",i,"\\]","\t","\\[",k,"\\]",sep=""),command)
              if (length(p) < 1){
                if (i != j){
                  a <- paste("[",i,"]","\t","[",j,"]","\t","[",k,"]","\t","[","gilv_GMS898_MX,gilv_JK03429_CA,plumb_GAV1335_MX","]","\t","##",sep="")
                  command <- append(command,a)
                }
              }
          }
        } 
      }  
    }   
}
setwd("~/Dropbox/vRAD3/phylo.analyses/Dtest/")
writeLines(command,"Dtest_allOuts.txt")

##summarizing output: 
#trim sample names to taxa names & paste together to make unique test id's (test = unique combination of taxa in p1/2+p3) 
df <- read.delim("output/out.D4.txt")
df$p1sp <- as.factor(gsub("\\[","",sapply(strsplit(as.character(df$P1),"_"),"[[",1)))
df$p2sp <- as.factor(gsub("\\[","",sapply(strsplit(as.character(df$P2),"_"),"[[",1)))
df$p3sp <- as.factor(gsub("\\[","",sapply(strsplit(as.character(df$P3),"_"),"[[",1)))
df$Osp <- as.factor(gsub("\\[","",sapply(strsplit(as.character(df$O),"_"),"[[",1)))
df$test <- as.factor(paste(df$p1sp,df$p2sp,df$p3sp,df$Osp))

#empty table to record summary stats
table <- data.frame(test=character(),n.loci=numeric(),p.disc=numeric(),D.range=character(),Z.range=character(),
                    sig.permut=character())

#loop over tests. Get min/max D & Z, average number of loci, average percent discordant loci, & number of significant tests.
for (i in levels(df$test)){
  a <- subset(df,test==i)
  test <- as.character(a$test[1])
  Dmin <- min(abs(a$D))
  Dmax <- max(abs(a$D))
  D.range <- paste("(",Dmin,",",Dmax,")",sep="")
  Zmin <- min(a$Z)
  Zmax <- max(a$Z)
  Z.range <- paste("(",Zmin,",",Zmax,")",sep="")
  n.permut <- nrow(a)
  #convert Z to p value, apply holm-bonferroni correction (n comparisons = n sample combinations per test)
  a$p.uncorrected <- 1-pnorm(a$Z)
  a$p.corrected <- p.adjust(a$p.uncorrected,method="holm")
  #find number of significant sample permutations out of total permutations. Note numeric below gives one-tailed p value cutoff. 
  sig.permut <- nrow(a[which(a$p.corrected < 0.005),])
  sig.fraction <- paste(sig.permut,"/",n.permut,sep="")
  n.loci <- mean(a$nloci)
  p.disc <- mean(a$pdisc)
  b <- data.frame(test,n.loci,p.disc,D.range,Z.range,sig.fraction)
  table <- rbind(table,b)
}

##output/write to file
table <- cbind(colsplit(table$test," ",names=c("P1","P2","P3","P4")),table[2:6])
#write.table(table,"D_out.txt",sep="\t")