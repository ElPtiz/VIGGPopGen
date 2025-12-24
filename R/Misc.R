require(rentrez)
require(rbibutils)


#' Technical internal func
#'
#' @param x
#'tbd
#' @returns
#'tbd
#' @examples
#'not needed
parse_journal<-function(x){
  pmid<-x$uid
  authors<-x$authors$name
  authors<-gsub(" ", ", ", authors)
  authors<-paste(authors, collapse=". and ")
  title<-x$title
  journal<-x$fulljournalname
  volume<-x$volume
  number<-x$issue
  pages<-x$pages
  doi<-strsplit(x$elocationid, " ")[[1]][2]
  rtype<-x$pubtype[1]
  year<-strsplit(x$pubdate, " ")[[1]][1]
  #if("Has Abstract"%in%x$attributes){
   #abstract<-entrez_fetch(db="pubmed", id=pmid, rettype="abstract")
  #}else{abstract<-NA}
  bib<-paste("@article{",pmid,",\nauthor = {", authors, "},\n",
             "title = {",title,"},\n",
             "year = {", year,"},\n",
             "journal = {",journal, "},\n",
             "volume = {", volume, "},\n",
             "number = {", number, "},\n",
             "pages = {", pages, "},\n",
             "doi = {", doi, "},\n",
             #"abstract ={", abstract, "},\n",
             "type = {", rtype, "},\n}", sep="")
  return(bib)
}

#' Technical internal func
#'
#' @param x
#'tbd
#' @returns
#'tbd
#' @examples
#'not needed
parse_summary_list<-function(x){
  reffile<-NA
  i<-1
  while(i<=length(x)){
    curr<-x[[i]]
    ref<-parse_journal(curr)
    reffile<-paste(reffile, ref, sep="\n")
    i=i+1
  }
  return(reffile)
}

#' Technical internal func
#'
#' @param x
#'tbd
#' @returns
#'tbd
#' @examples
#'not needed
parse_pmid_list<-function(x){
  reflist<-NA
  i1<-1
  i2<-200
  while(i1<=length(x)){
    currsum<-entrez_summary(db="pubmed", id=x[i1:i2])
    currparsed<-parse_summary_list(currsum)
    reflist<-paste(reflist, currparsed, sep="\n")
    message(c("Parsed ", i2, " articles out of ", length(x)))
    i1=i1+200
    i2=i2+200
  }
  write(reflist, "fetched_references.bib")
}

#' Search PubMed and export results as a .bib file
#'
#' Function utilizes Rentrez functional to search NCBI Pubmed and custom parser to export results into a .bib file ready to umport into Endnote or other soft.
#'
#' @param term
#' Search term as defined by Entrez search usage. Function can't handle too broad terms which return more then 9999 results.
#'
#' @returns
#' A .bib file in the working directory with formatted search results.
#' @export
#'
#' @examples
#' \dontrun{LitSearch("((animal migration[MeSH Terms] OR migrations[MeSH Terms]) AND(gene[All fields] OR gene[MeSH Terms] OR genetic[All fields] OR genetic[MeSH Terms] OR genetics[All fields] OR genetics[MeSH Terms]) AND (aves[ORGN] OR mammalia[ORGN]))")}
LitSearch<-function(term){

arts<- entrez_search(db="pubmed", term = term,
                     use_history = T, sort="pub_date")
artinf<-entrez_fetch(db="pubmed", web_history = arts$web_history, rettype = "uilist", parsed = F)
artinf<-gsub("[\r\n]", ",", artinf)
artinf<-strsplit(artinf, ",")
artinf<-artinf[[1]]
if(length(artinf)==9999){
  message("Retrieved 9999 articles, your search term is likely too broad, please refine")
  stop("Term too broad")
}
message(paste("Retrived information on ", length(artinf), " articles. Enter y to continue or n to stop and refine search term."))
cont<-readline("Continue? \n")
if(cont=="y"){
parse_pmid_list(artinf)}else{
  stop("Manually stopped")
}

}


