## loci2gphocs: converts pyRAD .loci output to sequence input file for gphocs
# default behavior is to output a list of character strings, one per line of the output file, use the write function to write to file (default settings work fine)
#example: 
#source("PathTo/loci2blat.R")
#output <- loci2gphocs("pathTo.lociFile")
#write(output,"pathToOutputFileLocation")

loci2gphocs <- function(infile){
  infile <- readLines(infile)
  a <- grep("//",loci)
  j <- 1
  output <- c(paste(length(a)),"")
  for (i in (1:length(a))){
    locus <- loci[j:(a[i]-1)]
    locus <- gsub(">","",locus)
    locus <- gsub("-","N",locus)
    j <- a[i]+1
    c <- strsplit(locus," [a-zA-Z]")
    d <- do.call(rbind.data.frame, c)
    nchar <- nchar(as.character(d[1,2]))
    nsamples <- nrow(d)
    header <- paste("locus",i," ",nsamples," ",(nchar+1),sep="")
    concat <- c(header,locus,"")
    output <- append(output,concat)
  }
  output
}

