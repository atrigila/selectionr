#' Subset species
#'
#' This function allows you to select only a few of the downloaded species to perform analysis.
#'
#' The function takes a table of species available in your MSA ranked by "priority" and selects those to create a new fasta file.
#'
#' @param gene.list Gene name (file name) of the MSA
#' @param orthologue.selection.directory Output directory
#' @param download.directory Original directory where files are stored
#' @param species_priority_table Table specifying which species should be prioritized
#' @keywords subset
#' @importFrom Biostrings readDNAStringSet
#' @importFrom utils write.table write.csv read.csv
#' @importFrom utils read.table
#' @export prioritize.species
#' @return A fasta file with the selected species.

prioritize.species <- function(gene.list, orthologue.selection.directory, download.directory,
                               species_priority_table = "reemplazos.txt") {

  setwd(orthologue.selection.directory)

  #order.species.to.keep <- "reemplazos.txt"
  #file.copy(paste0(software.directory, "/", order.species.to.keep), orthologue.selection.directory)

  especies.sirven <- data.frame(gene = character(),
                                spp.1 = character(),
                                spp.2 = character(),
                                spp.3 = character(),
                                spp.4= character(),
                                spp.5 = character(),
                                spp.6 = character(),
                                spp.7 = character(),
                                spp.8 = character(),
                                spp.9 = character(),
                                spp.10= character(),
                                spp.11 = character(),
                                spp.12 = character(),
                                spp.13 = character(),
                                spp.14 = character(),
                                stringsAsFactors =  FALSE)

  i <- 1


  priority.table <- read.table(species_priority_table, fill = TRUE, stringsAsFactors = FALSE)

  absent.orthologues <- data.frame(gene = character(), stringsAsFactors =  FALSE)
  enough.orthologues <- data.frame(gene = character(), stringsAsFactors =  FALSE)
  for (gene.name in gene.list) {

    setwd(download.directory)
    gen.file <- paste(gene.name, ".fasta", sep = "")
    dna <- readDNAStringSet(gen.file)
    tmp <- attributes(dna)$ranges
    species <- names(tmp)
    species


    sirven <- c()

    for (fila in c(1:13)) {
      for (columna in c(1:6)) {
        especie.tabla <- priority.table[fila,columna]
        if (especie.tabla %in% species) {
          sirven <- c(sirven, especie.tabla)
          break # para salir del ciclo
        }
      }
    }


    sirven.count <- length(sirven)
    rellenar.count <- 14-sirven.count
    relleno <- rep("", rellenar.count)
    fila <- c(gene.name, sirven, relleno)

    especies.sirven[i,] <- fila
    i <- i+1
  }

  setwd(orthologue.selection.directory)
  utils::write.csv(especies.sirven, file = "replacement_set.txt", row.names = FALSE) # this could be depleted and directly use species.sirven

  ## 3. Leer Fasta original, parsear secuencias solo por especies de reemplazo



  # priority.table <- read.table("reemplazos.txt", fill = TRUE, stringsAsFactors = FALSE)

 # priority.table # Tabla con especies de reemplazo

  spp.ok <- utils::read.csv("replacement_set.txt", stringsAsFactors = FALSE, na.strings = c("","NA")) # Leer especies presentes y faltantes por gen



  for (gen.index in 1:nrow(spp.ok)) {


    gen.1 <- spp.ok$gene[gen.index] # Seleccionar el gen 1 de la lista
    gen.2 <- paste(gen.1, ".fasta", sep = "") # Encontrar el archivo FASTA en el escritorio para ese gen

    fasta.read <-readDNAStringSet(paste0(download.directory, gen.2)) # Leer el archivo fasta completo de ese gen y cargarlo.
    salida <- data.frame(fasta = character(), stringsAsFactors = FALSE)

    for (i in 2:14) {
      current.spp <- spp.ok[gen.index,i]
      if (!is.na(current.spp))   {
        seq1 <- (fasta.read[current.spp])
        seq1.1 <- as.character(seq1)
        fasta <- paste(">", current.spp,"\n", seq1.1, sep="")
        salida[i-1,] <- c(fasta)
      }
    }
    fasta_human <- paste(">","homo_sapiens" ,"\n",as.character(fasta.read["homo_sapiens"]), sep="")
    salida <-  rbind(fasta_human,salida)

    if (nrow(salida) >= 14) {

      ##### Only write those genes with enough orthologues
      utils::write.table(salida, file = gen.2, quote = FALSE, col.names = FALSE, row.names = FALSE)
      last.enough <- nrow(enough.orthologues)
      new.enough <- last.enough+1
      enough.orthologues[new.enough,1]<-gen.1
      } else {
      last.error <- nrow(absent.orthologues)
      new.error <- last.error+1
      absent.orthologues[new.error,1]<-gen.1
      }
    }
  utils::write.table(absent.orthologues, file = "Not enough orthologues.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
  utils::write.table(enough.orthologues,  file = "Enough orthologues.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)

}
