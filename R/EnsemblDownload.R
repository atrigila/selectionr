#' Download orthologues from Ensembl
#'
#' An improved version to download orthologues from Ensembl.
#'
#' @param server Version of Ensembl to use
#' @param gene.name character; gene name to download (as symbol or Ensembl gene ID)
#' @param type character; one of "gene symbol" or "Ensembl"
#' @param target_species character; Use with gene symbol only.
#' @param target_taxon numeric; Taxon to download species.
#'
#' @return character stringr; A fasta file with sequences
#' @export download.orthologues
#'
#' @examples \dontrun{
#' x <- c("ENSAMXG00005000412", "ENSAMXG00005000517")
#' ex.genes <- lapply(setNames(x,x), download.orthologues,
#' server = "https://jan2020.rest.ensembl.org", type = "Ensembl",
#' target_taxon = 7898)}
#'
#'

download.orthologues <- function(server = "https://jan2020.rest.ensembl.org",
                                 gene.name,
                                 type = c("symbol", "Ensembl"),  target_species = NULL,
                                 target_taxon) {

  stopifnot(typeof(gene.name) == "character",
            typeof(target_taxon) == "double")

  type <- match.arg(type)
  if (type == "Ensembl") {
    if(!is.null(target_species)){
      warning(paste("Target species is not required when input is Ensembl gene id:", gene.name))
    }
    ext <- paste0("/homology/id/", gene.name,"?", 'target_taxon=',target_taxon, ";", 'sequence=cdna;type=orthologues')
  }
  if (type == "symbol") {
    stopifnot(typeof(target_species) == "character")
    ext <- paste0("/homology/symbol/",target_species,"/", gene.name, "?", 'sequence=cdna;type=orthologues')
  }

  response <- httr::GET(paste(server, ext, sep = ""), httr::content_type("application/json"))
  #httr::stop_for_status(response)

  if (is.null(response) | response$status_code == 400) {
    message(paste0("Null response from querying: ", gene.name))
    return("Empty")
  } else {
    table.of.query <- utils::head(jsonlite::fromJSON(jsonlite::toJSON(httr::content(response))))
    if (length(table.of.query$data$homologies[[1]]) == 0) {
      message(paste("No homologies for this gene:", gene.name))
      return("No homologies")
    } else {
      # output_name_FASTA_file <- paste0(download.directory,"/", gene.name, ".fasta")
      num.homologies <- length(table.of.query$data$homologies)
      for(num in num.homologies) {
        homologues <- table.of.query$data$homologies[[num]]
        condition.filter <- homologues$type=="ortholog_one2one"
        tabla.one2one <- homologues[condition.filter,]

        if(nrow(tabla.one2one) == 0) {
          message(paste("No 1-to-1 orthologs for this gene in this target taxon:", target_taxon, "for", num, "homology"))
          return("No homologies")
        } else {
          reference_species_table <- table.of.query$data$homologies[[num]]$source
          reference_species_name <- reference_species_table$species[1]
          reference_species_sequence <- reference_species_table$align_seq[1]
          reference_species_geneid <- reference_species_table$id[1]
          #  reference_species_protid <- reference_species_table$protein_id[1]
          #  reference_species_target_taxon_id <- reference_species_table$taxon_id[1]

          reference_species_as_1st_FASTA <- paste(">", reference_species_name,"\n", reference_species_sequence, sep="")
          result <- c(reference_species_as_1st_FASTA)

          tabla <- tabla.one2one$target
          columnas <- ncol(tabla)
          filas <-  nrow(tabla)
          for (i in 1:filas) {
            other_species_name <- tabla$species[i]
            other_species_sequences <- tabla$align_seq[i]
            #  other_species_geneid <- tabla$id[i]
            #  other_species_protid <-tabla$protein_id[i]
            #  other_species_taxonid <-tabla$taxon_id[i]

            each_species_fasta <- paste(">", other_species_name,"\n", other_species_sequences, sep="")
            result <- c(result, each_species_fasta)
            result <- gsub('-', '', result)
          }
        }
      }
      return(list(fasta = result))
    }
  }
  Sys.sleep(1)
}


