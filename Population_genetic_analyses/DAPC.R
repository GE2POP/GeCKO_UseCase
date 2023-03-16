# required packages 
library(vcfR)
library(adegenet)
library(ade4)
#library(pegas)

# read input data
args<-commandArgs(TRUE)
vcfFile=args[1] 
groupFile=args[2]
vcf<-read.vcfR(vcfFile)
# convert into genInd format
vcf_gi<-vcfR2genlight(vcf)
# add group/population information
group=read.table(groupFile, header=F,stringsAsFactors=TRUE)
colTitle<-c("group","genotype")
colnames(group)<-colTitle
group <- group[order(group$genotype),]
group$genotype
group$group

vcf_gi$pop<- group$group
vcf_gi$pop
#  K-means (with BIC criteria)
clust <- find.clusters(vcf_gi, max.n=10, scale=F, n.pca=118, parallel=F) 

# BIC plot
pdf(file="DAPC_BIC_plot.pdf")
plot(clust$Kstat, type="b", main="BIC", xlab="number of clusters (K)", ylab="BIC")
dev.off() 


# DAPC using K-means result
K <- clust$grp 

# DAPC with cross-validation
xval <- xvalDapc(vcf_gi, pop(vcf_gi), n.pca.max=118, n.rep=100, parallel = "multicore", ncpus = 5)
DAPC <- dapc(vcf_gi, pop=K, n.pca=10, n.da=5, parallel=F)  

# DAPC Plot
pdf("DAPC_plot.pdf")
FD1 <- DAPC$ind.coord[,1]
FD2 <- DAPC$ind.coord[,2]

plot(FD1,FD2, pch=16,cex=1,col=c("deepskyblue2","olivedrab3", "darkorange","violetred4")[pop(vcf_gi)], xlab="Dicriminant axis 1", ylab=" Dicriminant axis 2")
legend(30,-5 , legend=levels(pop(vcf_gi)), pch=16,col=c("deepskyblue2","olivedrab3", "darkorange","violetred4")) 
dev.off() 