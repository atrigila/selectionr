#' Species Pruner
#'
#' This function allows you to prune the general species tree according to the
#' species present in your alignment. It also generates a phylip-like file for
#' codeml input.
#'
#'
#' @param gene.list A list of multiple sequence alignments names
#' @param suffix Suffix/filetype of the files containing the multiple sequence
#'   alignments
#' @param input.directory Where alignments are stored
#' @param tree.output.directory Where trees are going to be written
#' @keywords trees
#' @importFrom phytools read.newick
#' @importFrom Biostrings readDNAStringSet
#' @importFrom ape drop.tip
#' @importFrom ape root
#' @importFrom ape unroot
#' @importFrom ape write.tree
#' @importFrom stats na.omit
#' @importFrom utils read.table write.table
#' @export create.custom.species.tree
#' @return A txt file containing trees and a phylip-like fasta file.



create.custom.species.tree <- function(gene.list, suffix = "ordered_nuc_trimmed.fas",
                                       input.directory, tree.output.directory){


  # Be careful here! Read the complete tree from the Ensembl version you are using


  setwd(tree.output.directory)
  tree <- phytools::read.newick(url("https://raw.githubusercontent.com/Ensembl/ensembl-compara/release/97/scripts/pipeline/species_tree.vertebrates.branch_len.nw"))
  tree$edge.length <- NULL

  for (gene.name in gene.list) {
    print(gene.name)
    setwd(input.directory)
    summary.ok <- data.frame("gene.name" = character(), stringsAsFactors = FALSE)
    summary.error <-  data.frame("gene.name" = character(), stringsAsFactors = FALSE)


    gene.name.omm <- paste0(input.directory, gene.name, suffix)
    dna <- Biostrings::readDNAStringSet(gene.name.omm)
    tmp <- attributes(dna)$ranges
    species <- names(tmp)
    species

    pruned.tree<-ape::drop.tip(tree, tree$tip.label[-stats::na.omit(match(species, tree$tip.label))])
    #ape::write.tree(pruned.tree)
    check <- length(species) == length(pruned.tree$tip.label)     #Sanity check species numb = tree tips numb = true
    check
    unrooted.pruned.tree <- ape::unroot(pruned.tree)
    destination.nwk <- paste(gene.name, '_pruned.tre', sep = "")
    file.destination.name <- paste0(tree.output.directory, destination.nwk)

    ape::write.tree(unrooted.pruned.tree, file = file.destination.name)

    setwd(input.directory)

    width.seq <- dna@ranges@width[1]
    length.seq <- length(species)
    width.seq
    length.seq

    if (width.seq %% 3 == 0) { # width.seq %% 3 means triplets ok
      table.gene.name <- utils::read.table(paste(gene.name.omm, sep = ""))
      max.index <- nrow(table.gene.name)+1
      temp <- data.frame(matrix(ncol=1, nrow=max.index))
      temp[1,1] <- as.character(paste("", length.seq, width.seq, sep = "\t"))
      temp[c(2:max.index),1] <- as.character(table.gene.name[,1])
      output.name <- paste(gene.name, ".phy", sep = "")
      phylip.file <- paste0(tree.output.directory,output.name)
      utils::write.table(temp, phylip.file, col.names = FALSE, row.names = FALSE, quote = FALSE, eol = "\n")

      print ("OK")
      print(output.name)

      summary.gene.ok <-  data.frame("gene.name" = as.character(gene.name),stringsAsFactors = FALSE)
      summary.ok <- rbind(summary.ok,summary.gene.ok)


    } else {
      warning(paste0("Total alignment length not %% 3 for: ", output.name))
      summary.gene.error <-  data.frame("gene.name" = as.character(gene.name),stringsAsFactors = FALSE)
      summary.error <- rbind(summary.error,summary.gene.error)
    }
  }




     datestamp <- date()
  write.table(summary.ok, file = paste0(datestamp," Genes_OK.txt"), row.names = FALSE)
  write.table(summary.error, file = paste0(datestamp," Genes_Error.txt"), row.names = FALSE)
}
