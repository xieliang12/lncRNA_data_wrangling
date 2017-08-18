#lncrnadb.org XML downloaded dataset handling
library(XML)
library(RCurl)
lncrnadb <- xmlParse("F:/RIBOBIO/LncRNA database/lncrnadb_org/Homo%20sapiens.xml")

#download each XML dataset file for each entry name of lncrnadb.org
name_path <- "//Name"
lnc_names <- sapply(lncrnadb[name_path], xmlValue)
url <- "http://www.lncrna.com/rest/search/"
replace_space_with_symbol <- function(x) x <- gsub(" ","%20",x)
replace_space_with_line <- function(x) x <- gsub(" ", "_", x)
lnc_names_urls <- unlist(lapply(lnc_names, replace_space_with_symbol))
files <- paste(lapply(lnc_names, replace_space_with_line), ".XML", sep="")
urls <- paste(url, lnc_names_urls, "/", sep="")
urls[106] <- "http://www.lncrna.com/rest/search/H19%20upstream%20conserved%201%20/"
files[8] <- "Alpha_250_Alpha_280.XML"
setwd("F:/RIBOBIO/LncRNA database/lncrnadb_org/")
#the code below is for download XML files using R function download.file
for (i in 1:length(urls)) {
  Sys.sleep(5)
  try(download.file(urls[i], files[i], mode="w"), silent=T)
}

get_annotation <- function(filename) {
  tmp <- xmlParse(filename)
  name_path <- "//Name"
  alias_path <- "//Alias"
  biotype_path <- "//Biotype"
  name <- sapply(tmp[name_path], xmlValue)
  alias <- sapply(tmp[alias_path], xmlValue)
  biotype <- sapply(tmp[biotype_path], xmlValue)
  annotation_name <- "//annotation//annotation[@section='Name']"
  annotation_gene <- "//annotation//annotation[@section='Gene copies']"
  annotation_character <- "//annotation//annotation[@section='Characteristics']"
  annotation_expression <- "//annotation//annotation[@section='Expression']"
  annotation_function <- "//annotation//annotation[@section='Function']"
  annotation_conservation <- "//annotation//annotation[@section='Conservation']"
  annotation_misc <- "//annotation//annotation[@section='Misc']"
  anno_name <- sapply(tmp[annotation_name], xmlValue)
  anno_gene_copies <- sapply(tmp[annotation_gene], xmlValue)
  anno_character <- sapply(tmp[annotation_character], xmlValue)
  anno_expression <- sapply(tmp[annotation_expression], xmlValue)
  anno_function <- sapply(tmp[annotation_function], xmlValue)
  anno_conservation <- sapply(tmp[annotation_conservation], xmlValue)
  anno_misc <- sapply(tmp[annotation_misc], xmlValue)
  if(is.character(name) && name != "") name else name <- "N/A"
  if(is.character(alias) && alias != "") alias else alias <- "N/A"
  if(is.character(biotype) && biotype != "") biotype else biotype <- "N/A"
  if(is.character(anno_name) && anno_name != "") {anno_name <- gsub("\n\n", "", anno_name)} else anno_name <- "N/A"
  if(is.character(anno_gene_copies) && anno_gene_copies != "") {anno_gene_copies <- gsub("\n\n","",anno_gene_copies)} else anno_gene_copies <- "N/A"
  if(is.character(anno_character) && anno_character != "") {anno_character <- gsub("\n\n","",anno_character)} else anno_character <- "N/A"
  if(is.character(anno_expression) && anno_expression != "") {anno_expression <- gsub("\n\n","",anno_expression)} else anno_expression <- "N/A"
  if(is.character(anno_function) && anno_function != "") {anno_function <- gsub("\n\n","",anno_function)} else anno_function <- "N/A"
  if(is.character(anno_conservation) && anno_conservation != "") {anno_conservation <- gsub("\n\n", "", anno_conservation)} else anno_conservation <- "N/A"
  if(is.character(anno_misc) && anno_misc != "") {anno_misc <- gsub("\n\n", "", anno_misc)} else anno_misc <- "N/A"
  df <- data.frame(name, alias, biotype, anno_name, anno_gene_copies, anno_character, 
            anno_expression, anno_function, anno_conservation, anno_misc)
  return(df)
}

get_species <- function(filename) {
  tmp <- xmlParse(filename)
  species <- sapply(tmp["//species"], xmlValue)
  if(species != "") {
    position_path <- "//species//Entry"
    name_path = "//Name"
    chr <- sapply(tmp[position_path], xmlValue)
    name <- sapply(tmp[name_path], xmlValue)
    if(is.character(name) && name != "") name else name <- "N/A"
    lnc_name <- rep(name, length(chr))    
    species_path <- "//species//Entry/@Species"
    species <- c()
    position <- c()
    for (i in 1:length(chr)) {
      species[i] <- tmp[species_path][[i]][[1]] 
      position[i] <- chr[i]
    }
    df <- data.frame(lnc_name, species, position)
    return(df)
  }
}

get_publication <- function(filename) {
  tmp <- xmlParse(filename)
  name <- sapply(tmp["//Name"], xmlValue)
  author_path <- "//Publication//Author"
  pubmed_path <- "//Publication//Pubmed"
  title_path <- "//Publication//Title"
  year_path <- "//Publication//Year"
  author <- sapply(tmp[author_path], xmlValue)
  pubmed <- sapply(tmp[pubmed_path], xmlValue)
  title <- sapply(tmp[title_path], xmlValue)
  year <- sapply(tmp[year_path], xmlValue)
  if(is.character(author) && is.character(name)) {
    names <- rep(name, length(author))
    df <- data.frame(names, author, pubmed, title, year)
  } else df <- NULL
  return(df)
}

get_sequence <- function(filename) {
  tmp <- xmlParse(filename)
  name <- sapply(tmp["//Name"], xmlValue)
  if(is.character(name) && name != "") name else name <- "N/A"
  sequence_record <- "//sequence//SequenceRecord/@record"
  species_path <- "//sequence//Species"
  fasta_path <- "//sequence//FastaSequence"
  accession_path <- "//sequence//Accession/@id"
  sequence_record_species <- sapply(tmp[species_path], xmlValue)
  names <- rep(name, length(sequence_record_species))
  accession_id <- tmp[accession_path]
  records <- c()
  species <- c()
  fasta <- c()
  accession <- c()
  if(is.character(sequence_record_species)) {
    for (i in 1:length(sequence_record_species)) {
    records[i] <- tmp[sequence_record][[i]][[1]]
    species[i] <- sequence_record_species[i]
    fasta[i] <- sapply(tmp[fasta_path], xmlValue)[i]
    if(!is.null(accession_id[[i]])) {accession[i] <- tmp[accession_path][[i]][[1]]} else accession[i] <- "N/A"
  }
  df <- data.frame(names, records, species, accession, fasta)
  return(df)
  } 
}

xml_files <- list.files(path="F:/RIBOBIO/LncRNA database/lncrnadb_org/lncrnadb/")
path <- paste("F:/RIBOBIO/LncRNA database/lncrnadb_org/lncrnadb/", xml_files, sep="")
lncrnadb_annotation <- data.frame(name=character(), alias=character(), biotype=character(), anno_name=character(),
                         anno_gene_copies=character(), anno_character=character(), anno_expression=character(),
                          anno_function=character(), anno_conservation=character(), anno_misc=character())
for (i in 1:length(path)) {
  lncrnadb_annotation <- rbind(lncrnadb_annotation, get_annotation(path[i])) 
}
lncrnadb_annotation$associated_with_cancer <- ""
lncrnadb_annotation$anno_function <- as.character(lncrnadb_annotation$anno_function)
for (i in 1:nrow(lncrnadb_annotation)) {
  if(grepl("cancer", lncrnadb_annotation[i,]$anno_function)) {
    lncrnadb_annotation[i,]$associated_with_cancer <- "yes"
  }
}
write.table(lncrnadb_annotation, file="F:/RIBOBIO/LncRNA database/lncrnadb_org/lncrnadb/lncrnadb_annotation.tsv",
            quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

lncrnadb_species <- data.frame(lnc_name=character(), species=character(), position=character())
for (i in 1:length(path)) {
  lncrnadb_species <- rbind(lncrnadb_species, get_species(path[i]))
}
write.table(lncrnadb_species, file="F:/RIBOBIO/LncRNA database/lncrnadb_org/lncrnadb/lncrnadb_species.tsv",
            quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

lncrnadb_publication <- data.frame(name=character(), author=character(), title=character(), year=character())
for (i in 1:length(path)) {
  lncrnadb_publication <- rbind(lncrnadb_publication, get_publication(path[i]))
}
write.table(lncrnadb_publication, file="F:/RIBOBIO/LncRNA database/lncrnadb_org/lncrnadb/lncrnadb_publication.tsv",
            quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

lncrnadb_sequence <- data.frame(names=character(), records=character(), species=character(), accession=character(),
                                fasta=character())
for (i in 1:length(path)) {
  lncrnadb_sequence <- rbind(lncrnadb_sequence, get_sequence(path[i]))
}

lncrnadb_sequence$accession <- as.character(lncrnadb_sequence$accession)
for (i in 1:nrow(lncrnadb_sequence)) {
  ifelse(grepl("&gt", lncrnadb_sequence[i,]$accession),
    lncrnadb_sequence[i,]$accession <- str_sub(lncrnadb_sequence[i,]$accession, (str_locate(lncrnadb_sequence[i,]$accession, "&gt;")[1,2]+1), 
            (str_locate(lncrnadb_sequence[i,]$accession, "&lt;/a&gt;")[1,1]-1)), lncrnadb_sequence[i,]$accession)
}
write.table(lncrnadb_sequence, file="F:/RIBOBIO/LncRNA database/lncrnadb_org/lncrnadb/lncrnadb_sequence.tsv",
            quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

