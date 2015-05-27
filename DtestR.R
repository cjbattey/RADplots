### Write D test input files in pyRAD format ###
#need to add extra commands block at top - see old files

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
            for(l in gilv){
              if (i != j){
                a <- paste("[",i,"]","\t","[",j,"]","\t","[",k,"]","\t","[",l,"]","\t","##",sep="")
                command <- append(command,a)
              }
            }
          }
        } 
      }  
    }   
}
setwd("~/Dropbox/vRAD3/phylo.analyses/Dtest/")
writeLines(command,"Dtest_all.txt")

