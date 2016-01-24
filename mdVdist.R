#mdVdist tests for a correlation bw missing data and genetic distance using a mantel test
#run time is proportional to the number of characters in each sequence... wouldnt recommend w/seqs > 100kbp. 
#to do: speed up pairwise md calc's with Biostrings package - strsplit+gsub too slow for long seqs.
#(also foreach(?) other way to get the pairwise N's? shooting for 1mbp in <1s...)

mdVdist <- function(file,nperm=1000,verbose=T,plot=T){
    require(data.table);require(ape);require(ade4)
  a <- data.frame(fread(file,header=T)) #read in sequences as a table
  names(a) <- gsub("X","",names(a)) #extract header info
  nchar <- as.numeric(names(a)[2]) #extract number of characters in sequence
  colnames(a) <- c("sample","seq")
  dna <- read.dna(file) #redundant read in for ape
  dna.dist <- dist.dna(dna,pairwise.deletion=T) #calculate dna distance (default model is k80)
  
  #calculate % shared loci (shared/total in alignment) for all pairwise comparisons
  md.dist <- apply(combn(nrow(a),2),2,function(x) {  
    p <- c(a$seq[x[1]],a$seq[x[2]]) #select a pair
    q <- strsplit(p,"") #split out each nucleotide
    q <- lapply(q,function(e) gsub("N",1,e)) #replace N's w/1
    q <- lapply(q,function(e) gsub("[^1]",0,e)) #replace everything else w/0
    q <- lapply(q,function(e) as.numeric(e)) 
    r <- q[[1]]+q[[2]]
    miss.none <- length(r[which(r == 0)])
    miss.one <- length(r[which(r == 1)])
    miss.two <- length(r[which(r == 2)])
    if(verbose==T){
      print(x)
      print(miss.one/(miss.none+miss.one))
    } 
    c(miss.none/(length(r)),x)
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

       
