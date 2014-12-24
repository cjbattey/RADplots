############################################# loci2blat(infile,outfile,cutsite) ##################################################
## Function to pick one sequence per locus from pyRAD's .loci output and write to a fasta file, for use with BLAT/BLAST
## infile=path to input .loci file; outfile=path to write output file; cutsite=complete cut site for restriction enzyme to add
## on to the beginning of the sequence. 

##NOTE: Possible problem with some loci for which pyRAD does not chop off the last base of the cut site - unclear why this is, but
##it results in an extra nucleotide when adding back on the complete cut site and will produce a gap in reference alignments. 
##Working on solutions.

loci2blat <- function(infile,outfile,cutsite){
lociLines <- readLines(infile)
endMarks <- grep("//",lociLines)
loci<- c()
j <- 1
for (i in 1:length(endMarks)){
  locus <- lociLines[j:endMarks[i]]
  j <- endMarks[i]+1
  n <- as.integer(runif(1,1,length(locus)-1))
  locus <- locus[n]
  locus <- sub(" {2,}",paste("\n",overhang,sep=""),locus)
  locus <- sub(">",paste(">locus",i,"_",sep=""),locus)
  loci <- append(loci,locus)
}
write(loci,outfile)
}

