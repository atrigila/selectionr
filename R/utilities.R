#' Export info from list as fasta
#'
#' Utility to export the sequences as
#' a fasta file. A fasta file is created for each element in the list
#'
#' @param list list; list obtained from download.orthologues be used as input
#' @param output.directory character; directory to write a fasta file for each
#'  element in the list
#' @param no.sequences logical; if \code{TRUE}, returns a txt file with a list of genes
#'  without sequences (Empty or No homologies)
#'
#' @return A fasta file
#' @export
#'
#' @examples \dontrun{
#' export.fasta.out(ex.genes, "/mydirectory/", no.sequences = FALSE)}
#'
export.fasta.out <- function(list, output.directory, no.sequences = TRUE) {
  stopifnot(dir.exists(output.directory))
  no.info.genes <- c()
  for (i in 1:length(list)) {
    if (list[[i]] == "No homologies" | list[[i]] == "Empty") {
      no.info.genes <- c(no.info.genes, names(list[i]))
    } else {
      fasta <- list[[i]]$fasta
      name.element <- names(list[i])
      utils::write.table(fasta, file = paste0(output.directory,"/", name.element,".fasta"))
    }
  }
  if(no.sequences == TRUE) {
    writeLines(no.info.genes, con = paste0(output.directory,"/", "No info genes.txt"))
  }
}

