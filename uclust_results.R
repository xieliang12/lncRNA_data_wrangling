#deal with the cluster results file output by uclust
results_uc <- read.table("F:/RIBOBIO/LncRNA database/three sources of LncRNA sequences and cluster/all_lncrna_results2.uc",
                       sep="\t", skip=8, stringsAsFactors=FALSE)
colnames(results_uc) <- c("type","cluster_nr","seq_length","percent_identity","strand",
                          "query_start","seed_start","alignment","query_label","target_label")
ind <- duplicated(results_uc[, c("cluster_nr", "query_label", "target_label")])
results_uc <- results_uc[!ind,]
results_uc_with_target <- subset(results_uc, target_label!="*")
results_uc_only_query <- subset(results_uc, target_label =="*")

grep_for_symbol <- function(df, item) {
  rows <- which(
  rowSums(
    `dim<-`(grepl(item, as.matrix(df), fixed=TRUE), dim(df))
    ) > 0
  )
}

gencode21_rows <- grep_for_symbol(results_uc_only_query, "ENST")
noncode_rows <- grep_for_symbol(results_uc_only_query, "NONHSAT")

#NONCODEv4 has 50004 unique lncRNAs
noncode <- results_uc_only_query[noncode_rows,]
noncode_unique <- noncode[!(noncode$query_label %in% results_uc_with_target$query_label),]
noncode_unique <- noncode_unique[!(noncode_unique$query_label %in% results_uc_with_target$target_label),]
cat("NONCODEv4 has", nrow(noncode_unique), "unique lncRNAs")

#gencode21 has 15203 unique lncRNAs
gencode21 <- results_uc_only_query[gencode21_rows,]
gencode21_unique <- gencode21[!(gencode21$query_label %in% results_uc_with_target$query_label),]
gencode21_unique <- gencode21_unique[!(gencode21_unique$query_label %in% results_uc_with_target$target_label),]
cat("gencode21 has", nrow(gencode21_unique), "unique lncRNAs")

#lncipedia has 10455 unique lncRNAs
gencode21_noncode_rows <- sort(c(gencode21_rows,noncode_rows))
lncipedia <- results_uc_only_query[setdiff(seq(1:nrow(results_uc_only_query)), gencode21_noncode_rows),]
lncipedia_unique <- lncipedia[!(lncipedia$query_label %in% results_uc_with_target$query_label),]
lncipedia_unique <- lncipedia_unique[!(lncipedia_unique$query_label %in% results_uc_with_target$target_label),]
cat("lncipedia has", nrow(lncipedia), "unique lncRNAs")

#find overlapped elements between two data sources
noncode_rows_target <- grep_for_symbol(results_uc_with_target, "NONHSAT")
gencode21_rows_target <- grep_for_symbol(results_uc_with_target, "ENST")
lncipedia_rows_target <- grep_for_symbol(results_uc_with_target, "lnc-")
length(intersect(noncode_rows_target, gencode21_rows_target))
#7037 lncRNAS exists in both NONCODEv4 and gencode21
length(intersect(noncode_rows_target, lncipedia_rows_target))
#6372 lncRNAs exists in both NONCODEv4 and lncipedia
length(intersect(gencode21_rows_target, lncipedia_rows_target))
#2925 lncRNAs exists in both gencode21 and lncipedia

results_uc_target <- results_uc_with_target[,c(2,9,10)]
agg_results_uc_target <- aggregate(results_uc_target[-3], by=list(results_uc_target$target_label),c)
colnames(agg_results_uc_target)[1] <- "target_label"
#get the row number found for three sources of lncRNA dataset
all_exist_rows <- intersect(intersect(grep_for_symbol(agg_results_uc_target, "ENST"), grep_for_symbol(agg_results_uc_target, "lnc-")), 
                            grep_for_symbol(agg_results_uc_target, "NONHSAT"))
#there are 2050 lncRNAs exist in all three sources of lncRNA datasets
cat("there are", length(all_exist_rows), "lncRNAs exist in all three sources of lncRNA datasets")

#get the query_label from results_uc_with_target dataset
results_uc_tq <- results_uc_with_target[,c(2,9)]
noncode_rows_to_remove <- results_uc_tq[grep_for_symbol(results_uc_tq, "NONHSAT"),]
lncipedia_rows_to_remove <- results_uc_tq[grep_for_symbol(results_uc_tq, "lnc-"),]
gencode_rows_to_remove <- results_uc_tq[grep_for_symbol(results_uc_tq, "ENST"),]

#add ">" to the front of the query_label
#add_arrow <- function(df, column) {
#  df[,column] <- paste(">", df[,column], sep="")
#  return(df)
#}

#write the to-remove sequence ids to three remove files respectively
#there are 14144, 6388, 7456 sequences to remove from noncodev4, lncipedia and gencodev21
write.table(noncode_rows_to_remove[,2], file="F:/RIBOBIO/LncRNA database/three sources of LncRNA sequences and cluster/to_remove_sequence_ids/noncodev4_remove.txt",
            quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(lncipedia_rows_to_remove[,2], file="F:/RIBOBIO/LncRNA database/three sources of LncRNA sequences and cluster/to_remove_sequence_ids/lncipedia_remove.txt",
            quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(gencode_rows_to_remove[,2], file="F:/RIBOBIO/LncRNA database/three sources of LncRNA sequences and cluster/to_remove_sequence_ids/gencode_remove.txt",
            quote=FALSE, row.names=FALSE, col.names=FALSE)

head(add_symbol(noncode_rows_to_remove, query_label))
results_uc_multiple_query <- results_uc_multiple_query[,c(2,3,9)]
library(splitstackshape)
LS <- cSplit(results_uc_multiple_query, 3, "|", direction="long")
