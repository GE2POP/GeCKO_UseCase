# required packages 
library(vcfR)
library(adegenet)
library(ade4)
library(ggplot2)


# read input data
args<-commandArgs(TRUE)
vcfFile=args[1] 
groupFile=args[2]

# read input vcf and convert into genInd format
vcf<-read.vcfR(vcfFile)
vcf_gi<-vcfR2genlight(vcf)

# add group/population information
group<-read.table(groupFile, header=F,stringsAsFactors=TRUE)
colTitle<-c("group","genotype")
colnames(group)<-colTitle
group <- group[order(group$genotype),]
group$genotype
group$group
group$group <- as.factor(group$group)
vcf_gi$pop<- group$group
vcf_gi$pop

# PCA with the first 3 axis
pca3=dudi.pca(df = vcf_gi, scannf = FALSE, nf = 3)

#  inertia
inertia=inertia.dudi(pca3,col.inertia = TRUE)
inertia$tot.inertia

# df for each combination
df1_2 <- data.frame(
  Individu = rownames(pca3$li),
  PC1 = pca3$li[, 1],
  PC2 = pca3$li[, 2],
  Groups = factor(pop(vcf_gi))
)

df1_3 <- data.frame(
  Individu = rownames(pca3$li),
  PC1 = pca3$li[, 1],
  PC3 = pca3$li[, 3],
  Groups = factor(pop(vcf_gi))
)

df2_3 <- data.frame(
  Individu = rownames(pca3$li),
  PC2 = pca3$li[, 2],
  PC3 = pca3$li[, 3],
  Groups = factor(pop(vcf_gi))
)

# PCA Plot
pdf("PCA_axis1_axis2.pdf")
ggplot(df1_2, aes(x = PC1, y = PC2, color = Groups, label = Individu)) +
  geom_point() +
  labs(x = "Axis 1", y = "Axis 2") + 
  theme_minimal() +
  scale_color_manual(values = c("deepskyblue2", "olivedrab3", "darkorange", "violetred4"),
                     labels = levels(pop(vcf_gi)))+
theme(legend.title = element_text(face = "bold"))
dev.off() 

pdf("PCA_axis1_axis3.pdf")
ggplot(df1_3, aes(x = PC1, y = PC3, color = Groups, label = Individu)) +
  geom_point() +
  labs(x = "Axis 1", y = "Axis 3") + 
  theme_minimal() +
  scale_color_manual(values = c("deepskyblue2", "olivedrab3", "darkorange", "violetred4"),
                     labels = levels(pop(vcf_gi)))+
theme(legend.title = element_text(face = "bold"))
dev.off()

pdf("PCA_axis2_axis3.pdf")
ggplot(df2_3, aes(x = PC2, y = PC3, color = Groups, label = Individu)) +
  geom_point() +
  labs(x = "Axis 2", y = "Axis 3") + 
  theme_minimal() +
  scale_color_manual(values = c("deepskyblue2", "olivedrab3", "darkorange", "violetred4"),
                     labels = levels(pop(vcf_gi)))+
theme(legend.title = element_text(face = "bold"))
dev.off()