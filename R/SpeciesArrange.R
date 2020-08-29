#' Species at top of the alignment
#'
#' A function to make a species be a top of the alignment
#'
#' Useful when some type of aligners rearrange species order. Works for DNA.
#'
#' @param gene.list character; The name of the gene
#' @param suffix character; The filename/type of the fasta file
#' @param species.atop character; The name of the species you want atop of the alignment (reference)
#' @param input.directory character; The input directory where the fasta file is stored
#' @param output.directory character; The output directory where you want the new fasta
#' @importFrom Biostrings readDNAStringSet
#' @importFrom Biostrings writeXStringSet
#'
#' @return Returns a new fasta file with your desired species at top of the alignment
#' @export species.top
#'
species.top <-  function(gene.list, suffix = "_nuc_trimmed.fas", species.atop = "homo_sapiens", input.directory, output.directory) {

  for (gene.name in gene.list) {
    gene.file <- paste0(input.directory, gene.name, suffix)
    dna <- Biostrings::readDNAStringSet(gene.file)
    moveToTop <- function(x, n) c(x[n], x[-n]) # https://support.bioconductor.org/p/127125/
    index <- which(species.atop == names(dna))
    new.order <- moveToTop(dna, index)
    output.name <- paste0(gene.name,"_ordered" ,suffix)
    Biostrings::writeXStringSet(x= new.order, filepath =  paste0(output.directory, output.name))
  }
}
