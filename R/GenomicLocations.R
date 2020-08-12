#' Get genomic location from Ensembl Protein
#'
#' This function allows you to get the geneomic location from aminoacid sites obtained from Ensembl proteins.
#'
#'
#' @param server Ensembl version and server to query
#' @param ensembl.protein.ID A list of Ensembl protein IDs
#' @param original.number Aminoacid number you would like to get its genomic location
#' @keywords genomic location
#' @importFrom httr GET
#' @importFrom httr stop_for_status
#' @importFrom jsonlite fromJSON
#' @importFrom jsonlite toJSON
#' @export ensembl.orthologue.download
#' @return A data frame with gene name, chr, start and end of original aminoacidic position.
#' @examples
#' \dontrun{genomic.locator(ensembl.protein.ID = "ENSP00000430945", original.number = 1269)}


genomic.locator <- function(server = "https://jul2019.rest.ensembl.org", ensembl.protein.ID,original.number){
  #library(httr)
  #library(jsonlite)

  # server <- "https://jul2019.rest.ensembl.org"
  ext <- paste0("/map/translation/", ensembl.protein.ID,"/")
  sites <- paste0(original.number,"..",original.number,"?")
  r <- GET(paste(server, ext, sites, sep = ""), content_type("application/json"))

 # stop_for_status(r)

  Sys.sleep(1)

  str <- head(fromJSON(toJSON(content(r))))
  chr <- paste0("chr",str$mappings$seq_region_name)
  start <- str$mappings$start[[1]]
  end <- str$mappings$end[[1]]

  df <- data.frame("gene.name"  = (ensembl.protein.ID),
    "chr" = c(chr),
             "start" = c(start),
             "end" = c(end), stringsAsFactors=FALSE)

  return(df)

}
