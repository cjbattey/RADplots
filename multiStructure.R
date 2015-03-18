##Function to write shell scripts for executing parallel runs of the program STRUCTURE.
##File paths can be complete or from the working diretory you intend to use.
##structureCommand is the command to run structure - either the path to the
#structure.bin file (it's in the .app - right click and go to show contents on a Mac to find it), or just "structure" if you have it in /usr/bin
##Output will be nruns shell scripts saved to separate lines in file threadme.txt in your working directory.
##
##To run these shell commands in parallel using the desired number of threads, use threader.py from https://github.com/slager/threader

multiStructure <- function(outfile,structureCommand,mainparams,nruns){
  shell <- c()
    for (run in 1:nruns){
      command <- paste(structureCommand," -m ",mainparams," -o ",paste(outfile,"_run",run,sep="")," > ",paste(outfile,"_run",run,".log",sep=""),sep="")
      shell <- append(shell,command)
    }
  write(shell,"threadme.txt")
}
