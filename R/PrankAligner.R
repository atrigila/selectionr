#' A PRANK aligner wrapper
#'
#' A function to align sequences with default prank settings.
#'
#' @param gene.name character; List of gene names (file names) to align
#' @param input.directory character; Directory where the fasta files with sequences
#'  are located
#' @param output.file.directory character; Directory where the fasta files with
#'  MSA will be written
#' @param software.directory character; Directory where your prank executable is located
#'
#' @return An aligned fasta file
#' @export prank.align
#' @examples \dontrun{prank.align(sample.alignment.2, input.directory = "~/data/",
#' output.directory = "~/data/", software.directory = "/software")}

prank.align <- function (gene.name, input.directory, output.file.directory, software.directory) {
  input.sequences <- paste0(input.directory, "/", gene.name, '.fasta')
  stopifnot(file.exists(input.sequences))
  file.copy(paste0(software.directory, "/prank"), output.file.directory)

   setwd(output.file.directory)
   prank.align <- paste0("./prank ","-d=", "'", input.sequences,"'"," -o=" , gene.name , "_pranked" , " -translate", " -F")
   system(prank.align)

   file.remove(paste0(output.file.directory, "/prank"))
}
