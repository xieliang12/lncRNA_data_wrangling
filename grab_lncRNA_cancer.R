library("RCurl")
library("XML")
library("plyr")

query.term <- c("long non-coding RNA and cancer", "lncRNA and cancer", "lincRNA and cancer", "long noncoding RNA and cancer")

pubmed_search <- function(query.term) {
  m <- ""
  for (i in 1:length(query.term)) {
    res <- EUtilsSummary(query.term[i], type="esearch", db="pubmed")
    fetch <- EUtilsGet(res)
    results <- data.frame(PMID(fetch), ArticleTitle(fetch), AbstractText(fetch))
    m <- rbind(m, results)
  }
  return(m)
}
results <- pubmed_search(query.term)
results <- data.frame(lapply(results, as.character), stringsAsFactors=FALSE)
results <- results[-1,]
results <- results[!duplicated(results$PMID.fetch.),]
names(results) <- c("pmid","title","abstract")
write.table(results, file="F:/RIBOBIO/LncRNA database/LncRNA disease database/ncbi_pubmed/lncRNA_cancer_by20150303.tsv",
            quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t", na="NA")
