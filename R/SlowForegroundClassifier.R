#' Slow foreground classifier
#'
#'Slow foreground classifier checks every single species at the alignment and
#' evaluates if they belong to a certain taxon or not. For several files, it is more
#' convenient to use the fast classifier.
#'
#' @param gene.list character; Gene name or file name
#' @param suffix character; Suffix of the file name
#' @param Taxon.id numeric; Taxon to check whether foreground species belong
#' @param input.directory character; Directory where files are located
#' @param output.directory character; Directory where csv files will be written
#'
#' @importFrom Biostrings readDNAStringSet
#' @importFrom stringr word
#' @importFrom httr GET content_type stop_for_status content
#' @importFrom jsonlite fromJSON toJSON
#' @importFrom utils write.csv2
#'
#' @return A csv file with species and class (either foreground or background)
#' @export species.foreground.classifier
#'
species.foreground.classifier <- function(gene.list, suffix = ".fasta", Taxon.id = "Mammalia",
                                          input.directory, output.directory) {

  for (gene.name in gene.list) {
      gen.start <- gene.name
      gen.end <- "_final_align_NT.aln"
      gen.file <- paste0(input.directory, gene.name, suffix)
      dna <- Biostrings::readDNAStringSet(gen.file)
      tmp <- attributes(dna)$ranges
      species <- names(tmp)
      species
      print(gene.name)

    classification.spp <- data.frame("names" = character(), "class" = character())

    for (specie in species) {

      specie.1 <- stringr::word(string = specie, start = 1, end = 1, sep = "_")
      server <- "https://rest.ensembl.org/taxonomy/classification/"
      individual.query <- paste0(server,specie.1, "?")
      print(individual.query)

      Sys.sleep(1)
      r <- httr::GET(paste(individual.query, sep = ""), httr::content_type("application/json"))

      httr::stop_for_status(r)

      ensembl.forspp <- (jsonlite::fromJSON(jsonlite::toJSON(httr::content(r))))

      if(Taxon.id %in% ensembl.forspp$scientific_name == "TRUE"){
        target.taxon.id.species <- data.frame("names" = specie, "class" = "foreground")
        classification.spp <- rbind(target.taxon.id.species, classification.spp)
      } else {
        target.taxon.id.species <- data.frame("names" = specie, "class" = "background")
        classification.spp <- rbind(target.taxon.id.species, classification.spp)
        }
      output.name <- paste0(gene.name, "_classification.csv")
      setwd(output.directory)
      utils::write.csv2(classification.spp, output.name, row.names = FALSE)
    }
  }
}
