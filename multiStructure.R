##Quick script to write a shell file for batch runs of structure via slagerthreader.py
##adjust runs and k range as needed.

runs <- c(1:10)
krange <- c(2:8)

shell <- c()
    for (i in krange){
      for(j in runs){
        randseed <- as.integer(runif(1,min=1,max=1000000))
        command <- paste("structure"," -m ",paste("./mainparams","_k",i,sep="")," -D ",
                         paste(randseed)," -o ",paste("./output/k",i,"_run",j,sep="")," > ",
                         paste("./output/k",i,"_run",j,".log",sep=""),"\n",sep="")
        shell <- append(shell,command)
      }
    }

write(shell,"structure.sh")
    
