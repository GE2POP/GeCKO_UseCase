# required packages  
library(ggplot2)
# EggLib input files 
DD10_DC10=read.table("egglib_stats_10DD_10DC.txt", header=TRUE, sep="\t", na.strings = "None")
DC10_DP10=read.table("egglib_stats_10DC_10DP.txt", header=TRUE, sep="\t", na.strings = "None")
DP10_DE10=read.table("egglib_stats_10DP_10DE.txt", header=TRUE, sep="\t", na.strings = "None")
# column selection with Fst values by population pairwise
DD10_DC10_ext<-DD10_DC10[,c("chrom", "start","FstW_DCxDD")]
DC10_DP10_ext<-DC10_DP10[,c("chrom", "start","FstW_DCxDP")]
DP10_DE10_ext<-DP10_DE10[,c("chrom", "start","FstW_DExDP")]
# format
DD10_DC10_ext<-as.data.frame(DD10_DC10_ext)
DC10_DP10_ext<-as.data.frame(DC10_DP10_ext)
DP10_DE10_ext<-as.data.frame(DP10_DE10_ext)
# Merge 
DD_DC_DP=merge(DD10_DC10_ext,DC10_DP10_ext, by="start")
DD_DC_DP_DE=merge(DD_DC_DP,DP10_DE10_ext, by="start")
# column selection with Fst values by population pairwise
fst<-DD_DC_DP_DE[,c("chrom.x","start","FstW_DCxDD","FstW_DCxDP","FstW_DExDP")]
colnames(fst)<-c("chrom","start","FstW_DCxDD","FstW_DCxDP","FstW_DExDP")
# selection of contigs close to the target genes
chr4B=subset(fst,chrom=="chr4B")
Rht=subset(chr4B,start>27500001)
RHT=subset(Rht,start<32500000)
chr4B_Rht=RHT[,c("chrom","start")]
chr4B_Rht$gene="Rht-B1b"

chr5A=subset(fst,chrom=="chr5A")
q=subset(chr5A,start>652500001)
Q=subset(q,start<657500000)
chr5A_Q=Q[,c("chrom","start")]
chr5A_Q$gene="Q"
gene=rbind(chr4B_Rht,chr5A_Q)
# merge
fst_gene=fst
fst_gene=merge(fst,gene, by="start", all.x=TRUE)
fst_gene$gene=as.character(fst_gene$gene)
fst_gene[which(is.na(fst_gene$gene)),"gene"]<-"autre"
fst_gene$gene=as.factor(fst_gene$gene)
# format
ESSAI_gene<-fst_gene
colnames(ESSAI_gene)<-c("start_contigs","chromosome_zavitan","FST_DD_DC","FST_DC_DP","FST_DP_DE","chrom.y","contigs")
ESSAI_gene<-ESSAI_gene[,c("chromosome_zavitan","start_contigs","contigs","FST_DD_DC","FST_DC_DP","FST_DP_DE")]

#--------------
### 4B : gene Rht

chr4B=subset(ESSAI_gene,chromosome_zavitan=="chr4B")
str(chr4B)
chr4B$FST_DD_DC=as.numeric(as.character(chr4B$FST_DD_DC))
chr4B$FST_DC_DP=as.numeric(as.character(chr4B$FST_DC_DP))
chr4B$FST_DP_DE=as.numeric(as.character(chr4B$FST_DP_DE))
str(chr4B)

#DD_DC (threshold value calculation + plot)
DD_DC_chr4B_order=na.omit(chr4B$FST_DD_DC[order(chr4B$FST_DD_DC)])
max(DD_DC_chr4B_order)
min(DD_DC_chr4B_order)

lengthDD_DC_chr4B_order=length(!is.na(DD_DC_chr4B_order))

valeurDD_DC_chr4B_0.95=0.95*lengthDD_DC_chr4B_order
seuilDD_DC_chr4B_0.95=DD_DC_chr4B_order[valeurDD_DC_chr4B_0.95]

mean_DD_DC_chr4B=mean(chr4B$FST_DD_DC, na.rm=TRUE)

pdf("Fst_scan_4B_DD_DC.pdf")
ggplot(data=chr4B, aes(x=start_contigs, y=FST_DD_DC, group=contigs))+
  geom_point(aes(color = contigs, shape = factor(contigs), size = contigs))+
  scale_color_manual(values=c("magenta3","grey70"))+
  scale_shape_manual(values=c(17, 19))+
  scale_size_manual(values=c(2, 1))+
  theme_light()+
  geom_hline(yintercept=seuilDD_DC_chr4B_0.95,linetype="dashed", color = "red")+
  scale_x_continuous(name="physical position", limits=c(0,870000000))+
  scale_y_continuous(name="FST", limits=c(-0.1,1))+
  theme(axis.text.x= element_text(size=15))+
  theme(axis.text.y= element_text(size=15))+
  theme(axis.title.x = element_text(size=18))+
  theme(axis.title.y = element_text(size=18))
dev.off()

#DC_DP (threshold value calculation + plot)
DC_DP_chr4B_order=na.omit(chr4B$FST_DC_DP[order(chr4B$FST_DC_DP)])
max(DC_DP_chr4B_order)
min(DC_DP_chr4B_order)

lengthDC_DP_chr4B_order=length(!is.na(DC_DP_chr4B_order))

valeurDC_DP_chr4B_0.95=0.95*lengthDC_DP_chr4B_order
seuilDC_DP_chr4B_0.95=DC_DP_chr4B_order[valeurDC_DP_chr4B_0.95]

mean_DC_DP_chr4B=mean(chr4B$FST_DC_DP, na.rm=TRUE)

pdf("Fst_scan_4B_DC_DP.pdf")
ggplot(data=chr4B, aes(x=start_contigs, y=FST_DC_DP, group=contigs))+
  geom_point(aes(color = contigs, shape = factor(contigs), size = contigs))+
  scale_color_manual(values=c("magenta3","grey70"))+
  scale_shape_manual(values=c(17, 19))+
  scale_size_manual(values=c(2, 1))+
  theme_light()+
  geom_hline(yintercept=seuilDC_DP_chr4B_0.95,linetype="dashed", color = "red")+
  scale_x_continuous(name="physical position", limits=c(0,870000000))+
  scale_y_continuous(name="FST", limits=c(-0.1,1))+
  theme(axis.text.x= element_text(size=15))+
  theme(axis.text.y= element_text(size=15))+
  theme(axis.title.x = element_text(size=18))+
  theme(axis.title.y = element_text(size=18))
dev.off()

#DP_DE (threshold value calculation + plot)
DP_DE_chr4B_order=na.omit(chr4B$FST_DP_DE[order(chr4B$FST_DP_DE)])
max(DP_DE_chr4B_order)
min(DP_DE_chr4B_order)

lengthDP_DE_chr4B_order=length(!is.na(DP_DE_chr4B_order))

valeurDP_DE_chr4B_0.95=0.95*lengthDP_DE_chr4B_order
seuilDP_DE_chr4B_0.95=DP_DE_chr4B_order[valeurDP_DE_chr4B_0.95]

mean_DP_DE_chr4B=mean(chr4B$FST_DP_DE, na.rm=TRUE)

pdf("Fst_scan_4B_DP_DE.pdf")
ggplot(data=chr4B, aes(x=start_contigs, y=FST_DP_DE, group=contigs))+
  geom_point(aes(color = contigs, shape = factor(contigs), size = contigs))+
  scale_color_manual(values=c("magenta3","grey70"))+
  scale_shape_manual(values=c(17, 19))+
  scale_size_manual(values=c(2, 1))+
  theme_light()+
  geom_hline(yintercept=seuilDP_DE_chr4B_0.95,linetype="dashed", color = "red")+
  scale_x_continuous(name="physical position", limits=c(0,870000000))+
  scale_y_continuous(name="FST", limits=c(-0.1,1))+
  theme(axis.text.x= element_text(size=15))+
  theme(axis.text.y= element_text(size=15))+
  theme(axis.title.x = element_text(size=18))+
  theme(axis.title.y = element_text(size=18))
dev.off() 

#------------
### 5A : gene Q

chr5A=subset(ESSAI_gene,chromosome_zavitan=="chr5A")
str(chr5A)
chr5A$FST_DD_DC=as.numeric(as.character(chr5A$FST_DD_DC))
chr5A$FST_DC_DP=as.numeric(as.character(chr5A$FST_DC_DP))
chr5A$FST_DP_DE=as.numeric(as.character(chr5A$FST_DP_DE))
str(chr5A)
#DD_DC (threshold value calculation + plot)
DD_DC_chr5A_order=na.omit(chr5A$FST_DD_DC[order(chr5A$FST_DD_DC)])
max(DD_DC_chr5A_order)
min(DD_DC_chr5A_order)

lengthDD_DC_chr5A_order=length(!is.na(DD_DC_chr5A_order))

valeurDD_DC_chr5A_0.95=0.95*lengthDD_DC_chr5A_order
seuilDD_DC_chr5A_0.95=DD_DC_chr5A_order[valeurDD_DC_chr5A_0.95]

mean_DD_DC_chr5A=mean(chr5A$FST_DD_DC, na.rm=TRUE)

pdf("Fst_scan_5A_DD_DC.pdf")
ggplot(data=chr5A, aes(x=start_contigs, y=FST_DD_DC, group=contigs))+
  geom_point(aes(color = contigs, size=contigs))+ 
  scale_color_manual(values=c("blue1","grey70"))+
  scale_size_manual(values=c(2,1))+
  theme_light()+
  geom_hline(yintercept=seuilDD_DC_chr5A_0.95,linetype="dashed", color = "red")+
  scale_x_continuous(name="physical position", limits=c(0,870000000))+
  scale_y_continuous(name="FST", limits=c(-0.1,1))+
  theme(axis.text.x= element_text(size=15))+
  theme(axis.text.y= element_text(size=15))+
  theme(axis.title.x = element_text(size=18))+
  theme(axis.title.y = element_text(size=18))
dev.off()

#DC_DP (threshold value calculation + plot)
DC_DP_chr5A_order=na.omit(chr5A$FST_DC_DP[order(chr5A$FST_DC_DP)])
max(DC_DP_chr5A_order)
min(DC_DP_chr5A_order)

lengthDC_DP_chr5A_order=length(!is.na(DC_DP_chr5A_order))

valeurDC_DP_chr5A_0.95=0.95*lengthDC_DP_chr5A_order
seuilDC_DP_chr5A_0.95=DC_DP_chr5A_order[valeurDC_DP_chr5A_0.95]

mean_DC_DP_chr5A=mean(chr5A$FST_DC_DP, na.rm=TRUE)

pdf("Fst_scan_5A_DC_DP.pdf")
ggplot(data=chr5A, aes(x=start_contigs, y=FST_DC_DP, group=contigs))+
  geom_point(aes(color = contigs, size=contigs))+ 
  scale_color_manual(values=c("blue1","grey70"))+
  scale_size_manual(values=c(2,1))+
  theme_light()+
  geom_hline(yintercept=seuilDC_DP_chr5A_0.95,linetype="dashed", color = "red")+
  scale_x_continuous(name="physical position", limits=c(0,870000000))+
  scale_y_continuous(name="FST", limits=c(-0.1,1))+
  theme(axis.text.x= element_text(size=15))+
  theme(axis.text.y= element_text(size=15))+
  theme(axis.title.x = element_text(size=18))+
  theme(axis.title.y = element_text(size=18))
dev.off() 

#DP_DE (threshold value calculation + plot)
DP_DE_chr5A_order=na.omit(chr5A$FST_DP_DE[order(chr5A$FST_DP_DE)])
max(DP_DE_chr5A_order)
min(DP_DE_chr5A_order)

lengthDP_DE_chr5A_order=length(!is.na(DP_DE_chr5A_order))

valeurDP_DE_chr5A_0.95=0.95*lengthDP_DE_chr5A_order
seuilDP_DE_chr5A_0.95=DP_DE_chr5A_order[valeurDP_DE_chr5A_0.95]

mean_DP_DE_chr5A=mean(chr5A$FST_DP_DE, na.rm=TRUE)

pdf("Fst_scan_5A_DP_DE.pdf")
ggplot(data=chr5A, aes(x=start_contigs, y=FST_DP_DE, group=contigs))+
  geom_point(aes(color = contigs, size=contigs))+ 
  scale_color_manual(values=c("blue1","grey70"))+
  scale_size_manual(values=c(2,1))+
  theme_light()+
  geom_hline(yintercept=seuilDP_DE_chr5A_0.95,linetype="dashed", color = "red")+
  scale_x_continuous(name="physical position", limits=c(0,870000000))+
  scale_y_continuous(name="FST", limits=c(-0.1,1))+  
  theme(axis.text.x= element_text(size=15))+
  theme(axis.text.y= element_text(size=15))+
  theme(axis.title.x = element_text(size=18))+
  theme(axis.title.y = element_text(size=18))
dev.off() 