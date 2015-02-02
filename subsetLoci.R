#########################################################subsetLoci.R##############################################################
###############################EXPERIMENTAL/may be full of errors. See especially line 21 and please tell me how to fix that##

################ Script to align RAD loci against a reference genome (via BLAT), 
################ subset for loci matching a selected chromosome (or chromosomes)
################ and produce a new alignment of only those loci.
################ Filters for blat hits are in the block beginning "#get an index of the locus # for z-chromosome loci"
################ Need pyRAD outfiles including one .phy-formatted file and the .loci format. Currently set up to work with sbf1 cut site RAD data. 

library(ggplot2);library(plyr)
setwd("/Users/fencelizard/Dropbox/vRAD3/")

#Reading in table of the number and size of chromosomes in the zebra finch - for optional plots of loci per chromosome. Ignore this for subsetting alignments.
z<- read.table("/Users/fencelizard/Documents/vRAD2/blat/zfinch_genome/zFinch_chrData.txt",header=TRUE)
z$Size <- as.numeric(as.character(z$Size))
z$ID <- factor(z$ID,as.character(z$ID))

#read in .loci file, find all bw loci breaks ("//"), for each locus read from the previous break +1 to the next break, 
#randomly choose one sample's sequence, store that sequence, add back on the full cut site (CCTGCAGG),
#reformat as fasta, name the loci sequentially, and write them to a new file to use with BLAT.
#PROBLEMS: sometimes pyRAD doesn't cut the whole overhang and instead leaves the last "G" on the beginning of the 
#sequence. Not clear why - it's not just loci with indels. Result is that adding the full cut site results 
#in a mismatch or gap of at least 1, b/c the new seq has an extra "G" at the end of the cut site.
lociLines <- readLines("./pyRAD/outfiles/c2h8.loci")
endMarks <- grep("//",lociLines)
loci<- c()
j <- 1
for (i in 1:length(endMarks)){
  locus <- lociLines[j:endMarks[i]]
  j <- endMarks[i]+1
  n <- as.integer(runif(1,1,length(locus)-1))
  locus <- locus[n]
  locus <- sub(" {2,}","\nCCTGCAGG",locus)
  locus <- sub(">",paste(">locus",i,"_",sep=""),locus)
  loci <- append(loci,locus)
}
write(loci,"./refAlign/seq/c2h8loci.fasta")

#run blat on the loci file (REMEMBER TO CHECK FILE PATHS FOR CORRECT FILES!)
system("blat /Users/fencelizard/Documents/vRAD2/blat/zfinch_genome/taeGut2.2bit /Users/fencelizard/Dropbox/vRAD3/refAlign/seq/c19h8loci.fasta /Users/fencelizard/Dropbox/vRAD3/refAlign/blatOutput/c19h8loci_zfAlign.pslx -t=dna -q=dna -out=pslx")

#get an index of the locus # for z-chromosome loci
a <- read.delim("./refAlign/blatOutput/c2h8loci_zfAlign.pslx",skip=5)
colnames(a) <- c("match","mis-match","rep.match","N's","Q.gap.count","Q.gap.bases","Tgap.count","Tgap.bases","strand","Q.name","Q.size","Q.start","Q.end","Tname","T.size","Tstart","T.end", "block.count","blockSizes","qStarts","tStarts","q.seq","t.seq")
a <- subset(a,a$match >=40)
duplicates <- subset(a,duplicated(a$Q.name))
a <- subset(a, !duplicated(a$Q.name))
a <- subset(a,a$Tname=="chrZ")
index <- as.character(a$Q.name)
index <- sapply(strsplit(index,"_"),'[[',1)
index <- sapply(strsplit(index,"s"),'[[',2)
index <- as.numeric(index)

######bar graphs for counts of loci by chromosome:
#a <- subset(a,a$Tname == "chr1"|a$Tname == "chr1A"|a$Tname == "chr1B"|a$Tname == "chr2"|a$Tname == "chr3"|a$Tname == "chr4"|a$Tname == "chr4A"|a$Tname == "chr5"|a$Tname == "chr6"|a$Tname == "chr7"|a$Tname == "chr8"|a$Tname == "chr9"|a$Tname == "chr10"|a$Tname == "chr11"|a$Tname == "chr12"|a$Tname == "chr13"|a$Tname == "chr14"|a$Tname == "chr15"|a$Tname == "chr16"|a$Tname == "chr17"|a$Tname == "chr18"|a$Tname == "chr19"|a$Tname == "chr20"|a$Tname == "chr21"|a$Tname == "chr22"|a$Tname == "chr23"|a$Tname == "chr24"|a$Tname == "chr25"|a$Tname == "chr26"|a$Tname == "chr27"|a$Tname == "chr28"|a$Tname == "chrZ"|a$Tname == "chrM")
#a$Tname <- as.character(a$Tname)
#a$Tname <- substr(a$Tname,4,nchar(a$Tname))
#a$min.samples <- "2"
#a$Tname <- factor(a$Tname,levels=c("1","1A","1B","2","3","4","4A","5",'6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','M','Z'))
#Loci <- c()
#for (i in levels(a$Tname)) {
#  sub <- subset(a,a$Tname == i)
#  rows <- nrow(sub)
#  Loci <- append(Loci,rows)
#}
#df <- data.frame(z,Loci)
#df$loci.per.base <- df$Loci/df$Size
#quartz()
#ggplot()+theme_bw()+scale_fill_grey()+
#  theme(axis.text.x=element_text(size=8),panel.grid.major = element_blank())+
#  ylab("Aligned Loci")+xlab("Chromosome")+
#  geom_bar(data=df,aes(x=ID,y=loci.per.base),stat="identity")

#Create empty data frame for alignment, subset out the z-loci, using the index vector from above to call the indices of the endMarks. Starting is hard: just have to input that.
#grab each locus, trim off the last line, replace spaces with a tab, 
lociLines <- readLines("./pyRAD/outfiles/c2h8.loci")
endMarks <- grep("//",lociLines)

samples <- read.table("./pyRAD/outfiles/c2h8.unlinked_snps",skip=1,colClasses="character")[1]
alignment <- data.frame(samples,seq=factor("\t"))
colnames(alignment) <- c('sample','seq')
lociLength <- c()
## this returns a list of vectors alternating sample,seq for each of the z Loci
for (i in index){
  locus <- lociLines[endMarks[i-1]:endMarks[i]]
  locus <- locus[3:length(locus)-1]
  locus <- gsub("* [-a-zA-Z]","\t",locus)
  locus <- sub("*\t[-a-zA-Z]","\tCCTGCAGG",locus)
  locus <- strsplit(locus,"\t")
  locus <- do.call(rbind.data.frame, locus)
  colnames(locus) <- c('sample','seq')
  locus$seq <- as.character(locus$seq)
  locus$sample <- as.character(locus$sample)
  locus$sample <- gsub(">","",locus$sample)
  locus$sample <- gsub("\t","",locus$sample)
  locus$sample <- gsub(" ","",locus$sample)
  lociLength <- append(lociLength,max(nchar(locus$seq)))
  alignment <- merge(alignment,locus,by="sample",all.x=TRUE)
}

#outputs a df with sample names and sequences, one column per locus, single NA in for missing data for a locus
#loop for each column: 1. read the length of the sequence block and store it, 2. search and replace NA with N repeated as many times as the sequence is long

for (i in 3:ncol(alignment)){
    a <- alignment[i]
    a[is.na(a)] <- paste(rep("N",lociLength[i-2]),collapse="")
    alignment[i] <-a
}

write.table(alignment,"zLoci_c2h8.phy")




