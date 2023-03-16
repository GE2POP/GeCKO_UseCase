# input file
args<-commandArgs(TRUE)
statFile=args[1]
tab=read.table(statFile, header=TRUE, sep="\t", na.strings = "None")
# column selection with Pi values (by contigs)
tab_short<-tab[,c("chrom", "start","lseff","DD_Pi","DC_Pi","DP_Pi","DE_Pi")]
# format
tab_short<-as.data.frame(tab_short)
# delete the 2 contigs on 2A and 7A that have a lseff<100bp x4
tab_short=subset(tab_short, lseff>100)
# Calculation of Pi values by site
tab_site_short=tab_short
tab_site_short$DD_Pi_site=(tab_site_short$DD_Pi)/(tab_site_short$lseff)
tab_site_short$DC_Pi_site=(tab_site_short$DC_Pi)/(tab_site_short$lseff)
tab_site_short$DP_Pi_site=(tab_site_short$DP_Pi)/(tab_site_short$lseff)
tab_site_short$DE_Pi_site=(tab_site_short$DE_Pi)/(tab_site_short$lseff)

# Calculation DRI values
summary(tab_site_short)
