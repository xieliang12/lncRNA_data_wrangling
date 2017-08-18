#cancer data cleaning up, before loading ng.3192-S4 save the excel file
#ng.3192-S4 to a tab-delimited .txt file
GSEA <- read.csv("E:/downloads/ng.3192-S4.txt", sep ="\t", skip=1)
GSEA <- GSEA[, c(1,4,5,8,9)]

NES <- aggregate(NES ~ Transcript.Name + Tissue.Type + Enrichment.Direction, data = GSEA, mean)
NES_sd <- aggregate(NES ~  Transcript.Name + Tissue.Type + Enrichment.Direction, data = GSEA, sd)
colnames(NES_sd)[4] <- "NES_Sd"
p_value <- aggregate(FWER.p.value ~ Transcript.Name + Tissue.Type + Enrichment.Direction, data=GSEA, mean)
p_value_sd <- aggregate(FWER.p.value ~ Transcript.Name + Tissue.Type + Enrichment.Direction, data=GSEA, sd)
colnames(p_value_sd)[4] <- "FWER.p.value_sd"

#merge multiple datasets#
merge_multiple <- function (x, y, x1, y1) {
  z <- merge(x, y,by=intersect(names(x), names(y)), all=TRUE)
  z1 <- merge(x1, y1, by=intersect(names(x1), names(y1)), all=TRUE)
  return(merge(z, z1, by=intersect(names(z), names(z1)), all=TRUE))
}

GSEA_mean_pvalue <- merge_multiple(NES, NES_sd, p_value, p_value_sd)
colnames(GSEA_mean_pvalue)[1:2] <- c("Transcript.Name","Tissue")
write.table(GSEA_mean_pvalue, file="F:/RIBOBIO/LncRNA database/lncRNA data clean up/GSEA_cancer_associated.csv",
          sep="\t", row.names=FALSE,  quote=FALSE, na="NA", col.names=TRUE)

#reference genome(hg19/GRCh37)
general_exp <- read.csv("E:/downloads/ng.3192-S3.txt", sep="\t", skip=1)
fpkm_columns <- c()
for (i in 1:length(colnames(general_exp))) {
  if (grepl("fpkm_95", colnames(general_exp)[i])) {
    fpkm_columns <- c(fpkm_columns, colnames(general_exp)[i])
    }
}
general_exp_subset <- cbind(general_exp[,c(1:23)], general_exp[, c(fpkm_columns)])

#merge the GSEA_mean_value and general_exp_subset table by column names Transcript_name and Tissue
expression_GSEA <- merge(general_exp_subset, GSEA_mean_pvalue, by=intersect(names(GSEA_mean_pvalue),names(general_exp_subset)),
                         all.x=TRUE)
expression_GSEA_cancer_lineage <- expression_GSEA[which(expression_GSEA$Association.Type=="Cancer and Lineage Associaiton"),]
for (i in 1:(nrow(expression_GSEA_cancer_lineage))) {
  expression_GSEA_cancer_lineage$match[i] <- paste0(expression_GSEA_cancer_lineage$Transcript.Name[i], "_",
                                         expression_GSEA_cancer_lineage$Tissue[i], "_", expression_GSEA_cancer_lineage$Enrichment.Direction[i])
}
expression_GSEA_cancer <- expression_GSEA[which(expression_GSEA$Association.Type=="Cancer Association"),]
for (i in 1:(nrow(expression_GSEA_cancer))) {
  expression_GSEA_cancer$match[i] <- paste0(expression_GSEA_cancer$Transcript.Name[i], "_",
                                                    expression_GSEA_cancer$Tissue[i], "_", expression_GSEA_cancer$Enrichment.Direction[i])
}
#the line below is to get rid of cancer association records that also found in cancer and lineage association datasets
expression_GSEA_cancer_nomatch <- expression_GSEA_cancer[which(!(expression_GSEA_cancer$match %in% expression_GSEA_cancer_lineage$match)),]
#subset the dataset that association_type is neither cancer_association nor cancer_lineage association
expression_GSEA_other <- expression_GSEA[which(!(expression_GSEA$Association.Type %in%  c("Cancer Association", "Cancer and Lineage Associaiton"))),]
expression_GSEA_all <- rbind(expression_GSEA_cancer_lineage[,-48], expression_GSEA_cancer_nomatch[,-48], expression_GSEA_other)
expression_GSEA_all <- expression_GSEA_all[order(expression_GSEA_all$Transcript.Name, expression_GSEA_all$Tissue),]
#clean the colnames of general_exp_subset a little bit, replace the last dot at the end of
#with nothing, replace the single or double dot with '_'
colnames_clean <- function (df) {
  for (i in 1:length(colnames(df))) {
  colnames(df)[i] <- gsub("\\.{1,}", "_", gsub("\\.$", "", colnames(df)[i]))
 }
}
write.table(expression_GSEA_all, file="F:/RIBOBIO/LncRNA database/lncRNA data clean up/expression_GSEA_all_hg19.csv", sep="\t",
          row.names=FALSE, col.names=TRUE, quote=FALSE)


#this part handles the expression data of lncRNA from the paper titled Human cancer long non-coding RNA transcriptomes.
normal_tissue <- read.csv("F:/RIBOBIO/LncRNA database/LncRNA disease database/cancer expression data/Normal tissue LncRNA expression.txt",
                          sep="\t", stringsAsFactors=FALSE, skip=2)
cancer_tissue <- read.csv("F:/RIBOBIO/LncRNA database/LncRNA disease database/cancer expression data/Cancer tissue lncRNA expression.txt",
                          sep="\t", stringsAsFactors=FALSE, skip=2)
common_columns <- intersect(names(normal_tissue), names(cancer_tissue))
normal_tissue <- normal_tissue[, common_columns]
cancer_tissue <- cancer_tissue[, common_columns]
#factor_to_character <- function (x) as.character(x)
#normal_tissue[, c(1:4)] <- lapply(normal_tissue[,c(1:4)], factor_to_character)
#cancer_tissue[, c(1:4)] <- lapply(cancer_tissue[,c(1:4)], factor_to_character)
add_normal_symbol <- function (x) paste0("normal_", x)
colnames(normal_tissue)[7:19] <- lapply(colnames(normal_tissue)[7:19], add_normal_symbol)
normal_plus_cancer <- merge(cancer_tissue, normal_tissue, by=intersect(names(cancer_tissue), names(normal_tissue)))


normal_plus_cancer_des <- normal_plus_cancer[,c(1,2,33,4,5,6)]
normal_plus_cancer <- normal_plus_cancer[, c(1,2,33,7:32)]
library(reshape2)
normal_plus_cancer_narrow <- melt(normal_plus_cancer, id=1:3)
normal_plus_cancer_narrow$variable <- as.character(normal_plus_cancer_narrow$variable)

normal_plus_cancer_narrow$type <- "cancer"
for (i in 1:nrow(normal_plus_cancer_narrow)) {
    if(grepl("normal_", normal_plus_cancer_narrow[i,]$variable, ignore.case=TRUE)) {
    normal_plus_cancer_narrow[i,]$type <- "normal"
    normal_plus_cancer_narrow[i,]$variable <- gsub("normal_", "", normal_plus_cancer_narrow[i,]$variable)
  }
}

names(normal_plus_cancer_narrow)[4:5] <- c("tissue", "TPM")
write.table(normal_plus_cancer_narrow, file="F:/RIBOBIO/LncRNA database/LncRNA disease database/cancer expression data/normal_vs_cancer_tpm.csv",
            quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

#read in the long non-coding RNA expression matrix for recalculation of all lncRNA
normal_vs_cancer <- read.csv("F:/RIBOBIO/LncRNA database/LncRNA disease database/cancer expression data/lncrna_expression.txt", 
                             sep="\t", stringsAsFactors=FALSE, header=TRUE)

select_columns <- c()
for (i in 1:length(colnames(normal_vs_cancer))) {
  if(grepl("brain|breast|embryonic_stem_cell|esophagus|gallbladder|liver|lung|lymph_node|muscle|prostate|retina|stomatch|thyroid",
     colnames(normal_vs_cancer)[i], ignore.case=TRUE)) select_columns <- c(select_columns, colnames(normal_vs_cancer)[i])
}

normal_vs_cancer <- cbind(normal_vs_cancer[,c(1:8)], normal_vs_cancer[,select_columns])
normal_vs_cancer_melt <- melt(normal_vs_cancer, id=1:8, variable.name="tissue", value.name="TPM", factorsAsStrings = TRUE)
colnames(normal_vs_cancer_melt)[9:10] <- c("tissue_detail","TPM")
normal_vs_cancer_melt$tissue_detail <- as.character(normal_vs_cancer_melt$tissue_detail)

library(stringr)
normal_vs_cancer_melt$type <- ifelse(str_detect(normal_vs_cancer_melt$tissue_detail,
                                                  ignore.case("Normal")), "normal", "cancer")
tissue_name <- function(x) return(strsplit(x, "_", fixed=TRUE)[[1]][1])
normal_vs_cancer_melt$tissue <- tolower(unlist(lapply(normal_vs_cancer_melt$tissue_detail, tissue_name)))
normal_vs_cancer_melt$tissue <- gsub("retinal","retina", normal_vs_cancer_melt$tissue)
normal_vs_cancer_melt$category <- paste0(normal_vs_cancer_melt$tissue, "_", normal_vs_cancer_melt$type)

#save the dataset for later use
write.table(normal_vs_cancer_melt, file="F:/R_data/normal_vs_cancer_melt.tsv", row.names=FALSE,
            col.names=TRUE, quote=FALSE, sep="\t")

normal_vs_cancer_melt <- read.csv("F:/R_data/normal_vs_cancer_melt.tsv", header=TRUE, sep="\t")
normal_vs_cancer_melt <- normal_vs_cancer_melt[, c(-9)]
colnames(normal_vs_cancer_melt)[1:8] <- c("ensembl_gene_id","gene_name","description", "chr","start","end","strand","biotype")


library(dplyr)
normal_vs_cancer_melt <- normal_vs_cancer_melt[which(normal_vs_cancer_melt$gene_name!=""),]
by_category <- group_by(normal_vs_cancer_melt, gene_name, category)
(grouped_TPM <- summarise(by_category, mean(TPM)))

library(reshape)
tissue_tpm_mean <- cast(grouped_TPM, gene_name ~ category, value="mean(TPM)")
tissue_tpm_mean <- tissue_tpm_mean[,c(1:5,7:24,6)]
tissue_tpm_mean <- tissue_tpm_mean[,c(1, seq(2,22,2), seq(3,23,2), 24)]
tissue_tpm_mean <- subset(tissue_tpm_mean, embryonic_normal > 0)

#the code below is for cancer expression to normal expressio fold change calculation#
newcolumns <- unlist(lapply(colnames(tissue_tpm_mean)[2:12], tissue_name))
colnames(tissue_tpm_mean)[25:35] <- newcolumns
for (i in 1:nrow(tissue_tpm_mean)) {
  for (j in c(2:12)) {
    if(tissue_tpm_mean[i,(j+11)]>0) {
      tissue_tpm_mean[i, (j+23)]=tissue_tpm_mean[i, j]/tissue_tpm_mean[i,(j+11)]
    } else
      tissue_tpm_mean[i, (j+23)]=NA
  }
}

annotation <- normal_vs_cancer_melt[,c(1:8)]
annotation <- annotation[!duplicated(annotation),]
merged <- merge(annotation, tissue_tpm_mean, by="gene_name", all.y=TRUE)
merged <- merged[!duplicated(merged),]


#put the HGNC synbol from Descripton column into a separate column#
merged$hgnc <- ""
for (i in 1:nrow(merged)) {
  if(grepl("hgnc", merged[i,]$description, ignore.case=TRUE)) {
    merged[i,]$hgnc <- gsub("\\]$", "", (strsplit(merged[i,]$description, split="Acc:")[[1]][2]))
    merged[i,]$description <- gsub(" \\[.*\\]$", "", merged[i,]$description)
  } 
}

write.table(merged, "F:/RIBOBIO/LncRNA database/LncRNA disease database/cancer expression data/recalculation_normal_vs_cancer.tsv",
            quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

