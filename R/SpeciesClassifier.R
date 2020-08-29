#' Foreground classifier
#'
#' Extract all unique species present in all alignments and assign to
#' foreground and background sets if they belong to a certain taxon id or not
#'
#' @param gene.list character; Gene list (file name) of MSAs to check for species names
#' @param suffix character; Suffix of the files (such as: ".fasta")
#' @param input.directory character; Directory where the files for the MSA are located
#' @param Taxon.id character; Foreground taxon to check if the species belongs to
#' @param classification.output.directory character; Directory where the csv will be written
#'
#' @importFrom Biostrings readDNAStringSet
#' @importFrom stringr word
#' @importFrom httr GET
#' @importFrom httr content_type stop_for_status content
#' @importFrom jsonlite fromJSON toJSON
#' @importFrom utils write.table write.csv2
#'
#' @return A csv with species names and class (either foreground or background)
#' @export species.classifier

species.classifier <- function(gene.list, suffix ="_nuc_trimmed.fas" , input.directory, Taxon.id = "Mammalia", classification.output.directory) {

all.species <- c()
  for (gene.name in gene.list) {
    gene.dna <- paste0(input.directory, gene.name,suffix)
    dna <- Biostrings::readDNAStringSet(gene.dna)
    species.per.gene <- names(dna)
    all.species <- c(all.species,species.per.gene)
  }

  species <- unique(all.species)


# Classify all unique species if they are foreground or background
# Ask Ensembl if the species genus belongs to my target foreground taxon id
  classification.spp <- data.frame("names" = character(), "class" = character())

  for (specie in species)  {

    #specie<- "homo_sapiens"
    specie.1 <- stringr::word(string = specie, start = 1, end = 1, sep = "_")
    server <- "https://rest.ensembl.org/taxonomy/classification/"
    individual.query <- paste0(server,specie.1, "?")
    print(individual.query)

    Sys.sleep(1)
    r <- httr::GET(paste(individual.query, sep = ""), httr::content_type("application/json"))

    httr::stop_for_status(r)
    Sys.sleep(1)
    pepe <- (jsonlite::fromJSON(jsonlite::toJSON(httr::content(r))))


    #Taxon.id <- "Mammalia"
    if(Taxon.id %in% pepe$scientific_name == "TRUE"){
      target.taxon.id.species <- data.frame("names" = specie, "class" = "foreground")
      classification.spp <- rbind(target.taxon.id.species, classification.spp)
    } else {
      target.taxon.id.species <- data.frame("names" = specie, "class" = "background")
      classification.spp <- rbind(target.taxon.id.species, classification.spp)
      }
    }


for (gene.name in gene.list) {

  gene.dna <- paste0(input.directory, gene.name,suffix)
  dna <- Biostrings::readDNAStringSet(gene.dna)
  species.per.gene <- names(dna)

  classification.spp.out <- data.frame("names" = character(), "class" = character())

  foreground.species <- classification.spp[classification.spp$class == "foreground", "names"]


  for (specie in species) {
    if(specie %in% foreground.species == "TRUE"){
      target.taxon.id.species <- data.frame("names" = specie, "class" = "foreground")
      classification.spp.out <- rbind(target.taxon.id.species, classification.spp.out)
    } else {
      target.taxon.id.species <- data.frame("names" = specie, "class" = "background")
      classification.spp.out <- rbind(target.taxon.id.species, classification.spp.out)
      }
    output.name <- paste0(classification.output.directory, gene.name, "_classification.csv")
    utils::write.csv2(classification.spp.out, output.name, row.names = FALSE)
    }
  }
}
