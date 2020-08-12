#' Write a fasta file from Ensembl Orthologues
#'
#' This non-exported function creates a fasta file from Ensembl orthologues information.
#'
#'
#' @param JSON.object JSON object to be parsed
#' @param out.file Output file name
#' @param name.gen Gene name to be parsed
#' @param out.file2 Not used parameter
#' @keywords writeFasta
#' @importFrom jsonlite fromJSON
#' @return A txt file containing sequences and species names in fasta format.


write.FASTA <- function(JSON.object, out.file, name.gen, out.file2)  {
  tabla_human <- JSON.object$data$homologies[[1]]$source
  spp <- tabla_human$species[1]
  sequence <- tabla_human$align_seq[1]
  human_seq <- paste(">", spp,"\n", sequence, sep="")


  result <- c(human_seq)

  tmp <- JSON.object$data$homologies[[1]]
  condition.filter <- tmp$type=="ortholog_one2one"
  tabla.one2one <- tmp[condition.filter,]

  tabla <- tabla.one2one$target



  columnas <- ncol(tabla)
  filas <-  nrow(tabla)

  if(filas > 0) {
    for (i in 1:filas) {
      spp <- tabla$species[i]
      sequence <- tabla$align_seq[i]

      temp <- paste(">", spp,"\n", sequence, sep="")
      result <- c(result, temp)
      result <- gsub('-', '', result) }
  } else {
    tabla_human <- JSON.object$data$homologies[[2]]$source
    spp <- tabla_human$species[1]
    sequence <- tabla_human$align_seq[1]
    human_seq <- paste(">", spp,"\n", sequence, sep="")
    result <- c(human_seq)
    tmp <- JSON.object$data$homologies[[2]]
    condition.filter <- tmp$type=="ortholog_one2one"
    tabla.one2one <- tmp[condition.filter,]

    tabla <- tabla.one2one$target
    for (i in 1:filas) {
      spp <- tabla$species[i]
      sequence <- tabla$align_seq[i]

      temp <- paste(">", spp,"\n", sequence, sep="")
      result <- c(result, temp)
      result <- gsub('-', '', result) }
  }



  ##### Prueba output tabla con IDs & species
  table.ids <- data.frame(tabla$species,tabla$protein_id)
  table.ids$gene.name <- rep(gene.name, nrow(table.ids))
  table.ids <- table.ids[,c(3,1,2)]


  print(result)
  write.table(result, out.file, quote = FALSE, row.names = FALSE, col.names = FALSE)
  # write.table(table.ids, out.file2, quote = FALSE, row.names = FALSE, col.names = TRUE)


  return(table.ids)

  getwd()
}



#' Download Ensembl Orthologues
#'
#' This function allows you to download 1-to-1 orthologues from Ensembl, given a target taxon.
#'
#'
#' @param gene.list A list of gene symbol names
#' @param target_taxon The NCBI number of the species you'd like to download
#' @param download.directory The directory where you'd like the downloaded multiple sequences to be stored
#' @keywords download
#' @importFrom Biostrings readDNAStringSet
#' @importFrom jsonlite fromJSON
#' @export ensembl.orthologue.download
#' @return A txt file containing sequences and species names in fasta format.
#' @examples
#' \dontrun{ensembl.orthologue.download(gene.list = c("RBFOX1", "FOXP2"), target_taxon = "9798", download.directory = ".")}


#   library(jsonlite)

ensembl.orthologue.download <-
function(gene.list, target_taxon, download.directory) {

  corrects <- data.frame(gene.name = character(), tabla.species = character(), tabla.protein_id = character (), stringsAsFactors = FALSE)
  incorrects <- data.frame(error.name = character(), stringsAsFactors = FALSE)



  print( ifelse( missing(gene.list), 'Gene list to download orthologues not specified', 'Correctly specified gene list' ) )
  print( ifelse( missing(target_taxon), 'Target taxon to search orthologues not specified', 'Correctly specified target taxon' ) )

  setwd(download.directory)




  for (gene.name in gene.list)
  {

    url.head <- 'https://jan2020.rest.ensembl.org/homology/symbol/human/'
    url.taxon <-paste('?target_taxon=',target_taxon,';', sep="")
    url.tail <- 'content-type=application/json;sequence=cdna;type=orthologues'
    destfile.custom <- paste(gene.name, '.json', sep="")
    url.custom <- paste(url.head, gene.name, url.taxon, url.tail, sep="")
    print(url.custom)

    JSON.object <- tryCatch(fromJSON(url.custom),
                            error = function (e) {
                              print(paste("This gene does not exist:", gene.name))
                              return(NULL)}
                              )
    salida <- paste(gene.name, ".fasta", sep ="")
    salida2 <- paste(gene.name, ".csv", sep ="")
    if (is.null(JSON.object)) {
      print ("JSON object is null")
      # HACER Error summary
      last.error <- nrow(incorrects)
      new.error <- last.error+1
      incorrects[new.error,1]<-gene.name

      } else {
      correct <- write.FASTA(JSON.object, salida, gene.name, salida2)
      corrects <- rbind(corrects,correct) }

    }
  Sys.sleep(1)

write.table(corrects, "DownloadedOrthologues.txt", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
write.table(incorrects, "ErrorDownloading.txt", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
}
