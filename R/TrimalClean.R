#' Trimal cleaning
#'
#'
#'Input aligned sequences as AA, trim them and then backtranslate to cds
#'
#' @param gene.list Gene name or file name
#' @param suffix Suffix for gene names
#' @param alignment.directory Directory where alignments are located
#' @param trim.output.directory Directory where trimmed sequences will be written
#' @param unaligned.directory Directory where unaligned sequences are located
#' @param software.directory Directory where Trimal executable is located
#'
#' @return A cleaned fasta file
#' @export trimal.clean

trimal.clean <- function(gene.list, suffix = "_pranked.best.fas", alignment.directory, trim.output.directory, unaligned.directory, software.directory) {

  file.copy(paste0(software.directory, "trimal"), alignment.directory)


  for (gene.name in gene.list) {
    setwd(alignment.directory)

    gene.alignment <- paste0(gene.name, suffix)
    out.name <- paste0(gene.name, "_aa_trimmed.fas")
    command.trim <- paste("./trimal","-in", gene.alignment,"-out", out.name, "-automated1" , sep = " ")
    system(command.trim)

    out.name.back <- paste0(gene.name, "_nuc_trimmed.fas")
    original.file <- paste0(unaligned.directory, "/", gene.name, ".fasta")

    command.back <- paste("./trimal","-in", out.name,"-out", out.name.back, "-backtrans", original.file , sep = " ")
    system(command.back)

    to.copy <-  list.files(path = alignment.directory, pattern = "_trimmed")

    file.copy(to.copy, trim.output.directory)
    file.remove(to.copy)

  }
  file.remove(paste0(alignment.directory,"trimal"))
}
