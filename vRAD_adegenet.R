setwd("~/Dropbox/vRAD3/phylo.analyses/adegenet/")
library(ggplot2);library(adegenet);library(reshape2)

###read in for different alignments (c[min#samplesPerLocus]h[maxSamplesWSharedHetSite].str)
#vireo <- read.structure("c4h8.str",col.lab=1,n.ind=35,n.loc=24098,col.pop=2,NA.char=-9,onerowperind=FALSE,ask=FALSE)
#vireo <- read.structure("c19h8.str",col.lab=1,n.ind=35,n.loc=7131,col.pop=2,NA.char=-9,onerowperind=FALSE,ask=FALSE)
#vireo<- read.structure("~/Dropbox/vRAD3/pyRAD/outfiles/",col.lab=1,n.ind=35,n.loc=7131,col.pop=2,NA.char=-9,onerowperind=FALSE,ask=FALSE)
#vireo <- read.structure("c38h8.str",col.lab=1,n.ind=35,n.loc=105,col.pop=2,NA.char=-9,onerowperind=FALSE,ask=FALSE)
vireo <- read.structure("c34h8.str",col.lab=1,n.ind=35,n.loc=938,col.pop=2,NA.char=-9,onerowperind=FALSE,ask=FALSE)
vireo <- read.structure("~/Dropbox/vRAD3/pyRAD/outfiles/c19h8.noOuts.str",col.lab=1,n.ind=35,n.loc=6201,col.pop=2,NA.char=-9,onerowperind=FALSE,ask=FALSE)
flav <- read.structure("flav/flav_c4h8.str",col.lab=1,n.ind=35,n.loc=938,col.pop=2,NA.char=-9,onerowperind=FALSE,ask=FALSE)

#pairwise fst
fst.vireo <- stamppConvert(vireo,"genlight")
pw.fst <- stamppFst(fst.vireo,nclusters = 4)

#pca
scaled.data <- scaleGen(vireo,missing="mean")
vireo.pca <- dudi.pca(scaled.data,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)
pdf("~/Dropbox/vRAD3/figures/c34h8_pca.pdf")
print(s.class(vireo.pca$li, pop(vireo),clabel=F,
        col=c("gold","forestgreen","magenta3","orangered","cornflowerblue","orange","sienna","dodgerblue4")))
dev.off()
#title("REVI PCA")
#add.scatter.eig(vireo.pca$eig[1:20], 3,1,2)


###dapc w/max #PC's, then test for opt # to keep w/alpha scores
temp <- dapc(vireo,n.pca=11,n.da=11)
opt.pc <- optim.a.score(temp)

###find clusters w/o using pop priors
clust <- find.clusters(vireo,max.n.clust=12)

find.clusters(vireo.small,max.n.clust=10)

###dapc and loadings
vireo@pop <- clust$grp
vireo@pop.names <- levels(clust$grp)

vireo.dapc <- dapc(vireo,n.pca=opt.pc$best,n.da=2)

###scatter and assignment plots
scatter(vireo.dapc,col=c("gold","chartreuse3","magenta2","deepskyblue","orangered"),scree.da=FALSE,clab=0)

assignplot(vireo.dapc)
compoplot(vireo.dapc,col=c("gold","chartreuse3","magenta","cornflowerblue","orangered"),lab="",
          txt.leg=c("altiloquus","flavoviridis","magister","olivaceus S","olivaceus N"))

###look at which SNP's differentiate pops
loadingplot(vireo.dapc$var.contr,threshold = 0.01)

#get list of names for strongly differentiating loci
disc.loci <- loadingplot(k$var.contr,threshold = 0.01)$var.names
disc.loci <- substr(disc.loci,1,4)
disc.loci <- disc.loci[which(duplicated(disc.loci)=="FALSE")]
#for each locus, calculate allele freq's in each population and store as a list of data frames
loci <- seploc(vireo)
loci.freqs <- list()
j <- 1
for (i in disc.loci){
  a <- tab(loci[[i]])
  b <- data.frame(apply(a, 2, function(e) tapply(e,pop(vireo),mean,na.rm=TRUE)))
  rownames(b) <- c("alt","flav","mag","olivN","olivS")
  loci.freqs[[j]] <- b
  j <- j+1
}
#plots for allele frequencies by population (change i for different loci)
i <- 8
df <- data.frame(pop=rownames(loci.freqs[[i]]),loci.freqs[[i]])
df <- melt(df,id.vars="pop")
ggplot(data=df,aes(x=pop,y=value,fill=variable))+theme_bw()+geom_bar(stat="identity")





