#' Get genomic location from Ensembl Protein
#'
#' This function allows you to get the geneomic location
#'  from aminoacid sites obtained from Ensembl proteins.
#'
#' @param server The REST Ensembl version you'd like to access
#' @param ensembl.protein.ID Ensembl Protein ID
#' @param original.number Original position to query
#'
#' @importFrom httr GET content_type content
#' @importFrom jsonlite fromJSON toJSON
#' @importFrom utils head
#'
#' @return A data frame with gene name, start and end in human genomic positions.
#' @export genomic.locator

genomic.locator <- function(server = "https://jul2019.rest.ensembl.org", ensembl.protein.ID, original.number){

  ext <- paste0("/map/translation/", ensembl.protein.ID,"/")
  sites <- paste0(original.number,"..",original.number,"?")
  r <- httr::GET(paste(server, ext, sites, sep = ""), httr::content_type("application/json"))
  Sys.sleep(1)

  str <- utils::head(jsonlite::fromJSON(jsonlite::toJSON(httr::content(r))))
  chr <- paste0("chr",str$mappings$seq_region_name)
  start <- str$mappings$start[[1]]
  end <- str$mappings$end[[1]]

  df <- data.frame("gene.name"  = (ensembl.protein.ID), "chr" = c(chr),
             "start" = c(start),"end" = c(end), stringsAsFactors=FALSE)
  return(df)
}
