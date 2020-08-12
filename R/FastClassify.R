#' Fast Foreground Classifier
#'
#' This function allows you to classify which of your sequences will be tested as foreground lineages.
#'
#' This function makes use of a list of species that you consider to be foreground and outputs a csv for each gene you are testing.
#' @param gene.list A list of multiple sequence alignments names
#' @param suffix Suffix/filetype of the files containing the multiple sequence alignments
#' @param foreground.species List of foreground species to be searched in the multiple sequences
#' @keywords foreground
#' @importFrom Biostrings readDNAStringSet
#' @export
#' @return A csv file containing species names and classification into foreground or background lineages. Useful for the \link[selectionr]{extract.tag.tree}.
#' @seealso \link[selectionr]{extract.tag.tree}
#' @examples
#' \dontrun{fast.classify(gene.list = c("RBFOX1", "FOXP2"), suffix = ".fasta", foreground.species = "homo_sapiens")}


## Fast foreground classifier

fast.classify <- function(gene.list,
                          suffix,
                          foreground.species)

for (gene.name in gene.list) {

  gen.start <- gene.name
  gen.file <- paste0(gene.name, suffix)
  dna <- Biostrings::readDNAStringSet(gen.file)
  tmp <- attributes(dna)$ranges
  species <- names(tmp)


  classification.spp <- data.frame("names" = character(), "class" = character())


  for (specie in species)
  {
    if(specie %in% foreground.species == "TRUE"){
      target.taxon.id.species <- data.frame("names" = specie, "class" = "foreground")
      classification.spp <- rbind(target.taxon.id.species, classification.spp)
    } else {
      target.taxon.id.species <- data.frame("names" = specie, "class" = "background")
      classification.spp <- rbind(target.taxon.id.species, classification.spp)

    }
    output.name <- paste0(gene.name, "_classification.csv")
    write.csv2(classification.spp, output.name, row.names = FALSE)

  }
}
