##the code below can handle the fasta sequence file with package seqinr##
install.packages("seqinr")
library(seqinr)
sequences <- read.fasta(file="F:/RIBOBIO/LncRNA database/LncRNA disease database/LncRNADisease/rna_seq_clean.txt", as.string=TRUE)
name <- c()
annotation <- c()
length <- c()
seq <- c()

for (i in 1:length(sequences)) {
  name[i] <- getName(sequences[[i]])
  annotation[i] <- getAnnot(sequences[[i]])
  length[i] <- getLength(sequences[[i]])
  seq[i] <- getSequence(sequences[[i]], as.string=TRUE)[[1]][1]
}

rna_seq <- data.frame(name, annotation, length, seq)

#extract the ncbi genebank id from name
extract_genbank <- function(x) return(str_split(x, "\\|")[[1]][4])
rna_seq$name <- unlist(lapply(rna_seq$name, extract_genbank))
####################################################
#load lncRNAdisease database file
lncrna_disease <- read.csv("F:/RIBOBIO/LncRNA database/LncRNA disease database/LncRNADisease/experimentally_supported LncRNA-disease.txt",
                           sep="\t", header=FALSE, stringsAsFactors=FALSE)
lncrna_disease <- lncrna_disease[,c(2,3,5:13)]
colnames(lncrna_disease) <- c("gene_name","disease","description","chr","start","end","strand","species","alias","genbank_id", "pmid")
lncrna_disease <- lncrna_disease[!duplicated(lncrna_disease),]

#output a bed file for sequences extract from hg19
bed <- function(df, sp) {
  df_name <- subset(df, start != "N/A" & species == sp)
  df_name$score <- "."
  df_name <- df_name[,c(4,5,6,1,12,7)]
  df_name <- df_name[!duplicated(df_name),]
}

human_lncrna_disease_bed <- bed(lncrna_disease, "Human")
mouse_lncrna_disease_bed <- bed(lncrna_disease, "Mus")

write.table(human_lncrna_disease_bed, "F:/RIBOBIO/LncRNA database/LncRNA disease database/LncRNADisease/human_lncrna_disease_bed",
            quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t") 
write.table(mouse_lncrna_disease_bed, "F:/RIBOBIO/LncRNA database/LncRNA disease database/LncRNADisease/mouse_lncrna_disease_bed",
            quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

lncrna_disease_seq <- read.csv("F:/RIBOBIO/LncRNA database/LncRNA disease database/LncRNADisease/lncrna_disease_seq.tab",
                               sep="\t", header=FALSE, stringsAsFactors=FALSE)
colnames(lncrna_disease_seq) <- c("position","seq")
human_lncrna_disease_bed
