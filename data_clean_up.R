#lncRNADisease database#
experi_support <- read.csv("F:/RIBOBIO/LncRNA database/LncRNA disease database/LncRNADisease/experimentally_supported LncRNA-disease.txt",
                           sep="\t", header=FALSE)
names(experi_support) <- c("id","lnc_name","disease","dysfunction","description","chr","start","end","strand",
                           "species","alias","genbank","pubmedid")

#load NONCODEv4u1_human_lncRNA and expression data
NONCODEv4u1 <- read.csv('F:/RIBOBIO/LncRNA database/LncRNA disease database/NONCODEv4/')

#load cluster lncRNA sequences#
results.uc <- read.csv("F:/Program Files (x86)/PuTTY/results.uc", skip=8, sep="\t", stringsAsFactors=FALSE)
names(results.uc) <- c("type","clusternum","Seqlength","identify","strand","QueryStart","SeedStart","Alignment",
                       "QueryLabel","TargetLabel")
results.uc$TargetLabel <- ifelse(results.uc$TargetLabel=="*", NA, results.uc$TargetLabel)
results.uc$QueryLabel <- ifelse(results.uc$QueryLabel=="*", NA, results.uc$QueryLabel)
results.uc$cluster <- NA
for (i in 1:nrow(results.uc)) {
  if((!is.na(results.uc[i,]$QueryLabel)) && (!is.na(results.uc[i,]$TargetLabel))) {
    results.uc[i,]$cluster <- 1 
  } 
}
results.uc.alignment <- subset(results.uc, cluster==1)
results.uc.alignment <- results.uc.alignment[order(results.uc.alignment$identify),]
nrow(results.uc.alignment)
