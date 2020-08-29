#' Prequal cleaner wrapper
#'
#' @param gene.list character; A list of gene names or file names.
#' @param suffix character; Suffix for the file names (".fasta", etc)
#' @param unaligned.directory character; Directory where unaligned sequences are stored
#' @param trim.output.directory character; Directory where cleaned sequences will be written
#' @param software.directory character; Directory where prequal is located
#'
#' @return A fasta file with trimmed sequences
#' @export prequal.clean

prequal.clean <- function(gene.list, suffix = ".fasta", unaligned.directory, trim.output.directory, software.directory) {
  file.copy(paste0(software.directory, "prequal"), trim.output.directory)
  for (gene.name in gene.list) {
    gene.sequences <- paste0(gene.name, suffix)
    file.copy(paste0(unaligned.directory, gene.sequences), trim.output.directory)
    setwd(trim.output.directory)

    command.trim <- paste("./prequal", "-noremoverepeat", gene.sequences, sep = " ")
    system(command.trim)
    file.remove(paste0(trim.output.directory, gene.sequences))
  }
  file.remove(paste0(trim.output.directory,"prequal"))
}
