##Function to write shell scripts for executing parallel runs of the program STRUCTURE. 
##File paths can be complete or from the working diretory you intend to use. 
##structureCommand is the command to run structure - either the path to the structure.bin file (it's in the .app - right click and go to show contents on a Mac to find it)
##Output will n shell scripts saved to your working directory, where n is the number of threads.
##
##To run these, just open a separate terminal window for each thread and run one shell script per terminal window. 
##Note that you'll also need to run <chmod +x shellscript.sh> to make the scripts executable before running them. 
##This is a pain in the butt, but I don't know a way around it. 

multiStructure <- function(outfile,structureCommand,mainparams,threads,nruns){ 
  shell <- c()
  run <- 1
  for (i in 1:threads) {
    for (j in 1:(nruns/threads)){
      command <- paste(structureCommand," -m ",mainparams," -o ",paste(outfile,"_run",run,sep="")," > ",paste(outfile,"_run",run,".log",sep=""),"\n",sep="")
      shell <- append(shell,command)
      run <- run+1
    }
    write(shell,paste("thread_",i,".sh",sep=""))
    shell <- c()
  }
}

