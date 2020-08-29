#' Tag trees
#'
#' This function allows you to tag trees in foreground or background lineages.
#'
#'
#' @param gene.list character; A list of multiple sequence alignments names
#' @param m0.directory character; The directory where your PAML m0 files are located
#' @param classifier.directory Tcharacter; he directory where csv files indicating which species are foreground and background.
#' @param Ha.directory character; The directory where the alternative hypotesis of your PAML will be run
#' @param Ho.directory character; The directory where the null hypotesis of your PAML will be run
#' @keywords tag.trees
#' @importFrom phytools read.newick
#' @importFrom Biostrings readDNAStringSet
#' @importFrom ape read.tree
#' @importFrom stringr str_extract
#' @importFrom stringr str_replace
#' @importFrom ape write.tree
#' @importFrom ape makeNodeLabel
#' @importFrom utils read.csv
#' @export extract.tag.tree
#' @return A txt file in the Ha and Ho directories containing the trees tagged with #1 in the foreground lineage.

extract.tag.tree <- function(gene.list, m0.directory, classifier.directory,
                             Ha.directory, Ho.directory) {

  setwd(m0.directory)
  for (gene.name in gene.list) {
    if (file.exists(paste0("mlc_M0_", gene.name)) == FALSE) {
      stop("File does not exist.")
    } else {  # Extract tree
    lines <- readLines(paste0("mlc_M0_", gene.name))
    lines.start.array <- startsWith(lines, "tree length =")
    pos.start <- which(lines.start.array == TRUE)
    lines.end.array <- startsWith(lines, "Detailed output identifying parameters")
    pos.end <- which(lines.end.array == TRUE)
    lineas.utiles <- lines[pos.start:pos.end]
    pos.enters <- which(lineas.utiles == "")
    utiles2 <- lineas.utiles[pos.enters[2]:pos.enters[3]]
    utiles2 <- paste(utiles2, sep = "", collapse = "")
    utiles2 <- trimws(utiles2)

    #Tag tree in foreground
    tr <- ape::read.tree(text = utiles2)
    #plot(tr, show.node.label = TRUE)
    tr3 <- tr

    all.spp <- tr$tip.label

    output.classified.gene <- paste0(classifier.directory,gene.name, "_classification.csv")
    csv <- utils::read.csv(output.classified.gene, sep = ";")
    colnames(csv) <- c("names", "class")
    csv[,"names"] <- tolower(csv[,"names"])
    csv <- csv[,c("names", "class")]

    fore.spp <- c()

    for  (species in all.spp) {
      if (csv[which(csv$names == species), 2] == "foreground") {
        fore.spp <-  c(fore.spp, species)
      }
    }

    if (length(fore.spp) == 1) {
      pattern <- paste0(fore.spp,": ", "[:digit:]", ".", "[:digit:]+")
      # pattern <- paste0("homo_sapiens: ", "[:digit:]", ".", "[:digit:]+")
      # str_view(utiles2, pattern)
      to.tag <- stringr::str_extract(utiles2, pattern)
      replacement <- paste0(to.tag, " #1")
      utiles3 <- stringr::str_replace(utiles2, to.tag, replacement)
      setwd(m0.directory)
      output.name <- paste0("tagged_mlc_M0_", gene.name, ".txt")
      writeLines(utiles3, con = output.name)
      file.copy(output.name, Ha.directory)
      file.copy(output.name, Ho.directory)
      file.copy(paste0(gene.name, ".phy"), Ha.directory)
      file.copy(paste0(gene.name, ".phy"), Ho.directory)
    } else {
      tr3 <- ape::makeNodeLabel(tr3, "u", nodeList = list("#1" = fore.spp))
     # plot(tr3, show.node.label = TRUE)
      setwd(m0.directory)
      output.name <- paste0("tagged_mlc_M0_", gene.name, ".txt")
      ape::write.tree(tr3, file = output.name)

      file.copy(output.name, Ha.directory)
      file.copy(output.name, Ho.directory)
      file.copy(paste0(gene.name, ".phy"), Ha.directory)
      file.copy(paste0(gene.name, ".phy"), Ho.directory)
      }
    }
  }
}
