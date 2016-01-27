#mdVdist tests for a correlation bw missing data and genetic distance using a mantel test
#run time is proportional to the number of characters in each sequence... wouldnt recommend w/seqs > 100kbp. 
#to do: speed up pairwise md calc's with Biostrings package - strsplit+gsub too slow for long seqs

mdVdist <- function(file,nperm=1000,verbose=T,plot=T){
  require(data.table);require(ape);require(ade4)
  a <- data.frame(fread(file,header=T)) #read in sequences as a table
  names(a) <- gsub("X","",names(a)) #extract header info
  nchar <- as.numeric(names(a)[2]) #extract number of characters in sequence
  colnames(a) <- c("sample","seq")
  dna <- read.dna(file) #redundant read in for ape
  dna.dist <- dist.dna(dna,pairwise.deletion=T) #calculate dna distance (default model is k80)
  
  #calculate % shared loci for all (unique) pairwise combinations of samples
  md.dist <- apply(combn(nrow(a),2),2,function(x) {  
    p <- c(a$seq[x[1]],a$seq[x[2]]) #select a pair
    p <- gsub("N",1,p) #replace N with 1
    p <- gsub("[^1]",0,p) #replace everything else with 0
    q <- strsplit(p,"") #split out each nucleotide
    q <- lapply(q,function(e) as.numeric(e)) #reformat to numeric
    r <- q[[1]]+q[[2]] #add the pair of converted sequences
    miss.none <- length(r[which(r == 0)])
    miss.one <- length(r[which(r == 1)])
    miss.two <- length(r[which(r == 2)])
    if(verbose==T){
      print(x)
      print(miss.one/nchar)
    } 
    c(miss.none/nchar,x)
  })

  #reformat list output to a matrix.
  m <- matrix(NA,nrow = nrow(a),ncol = nrow(a))
  for(i in 1:ncol(md.dist)){
    m[md.dist[,i][3],md.dist[,i][2]] <- md.dist[,i][1]
    m[md.dist[,i][2],md.dist[,i][3]] <- md.dist[,i][1]
    m[md.dist[,i][2],md.dist[,i][2]] <- 1
  }
  m <- as.dist(m)
  
if(plot==T){
  line <- lm(m~dna.dist)
  plot(x=dna.dist,y=m,ylab="% Shared Loci",xlab="Genetic Distance")+abline(line)
} 
  
  mantel.test(as.matrix(dna.dist),as.matrix(m),nperm=nperm)
}

       
