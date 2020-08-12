#!/usr/bin/env Rscript

# Species pruner according to species present in alignment
#install.packages("phytools")
#install.packages("ape")


create.custom.species.tree <- function(gene.list, suffix = "ordered_nuc_trimmed.fas", input.directory, tree.output.directory)
{
  library("phytools")
  library(Biostrings)
  library(ape)
  library(httr)
  library(stringr)
  
  # Be careful here! Read the complete tree from the Ensembl version you are using
  
  
  setwd(tree.output.directory)
  tree <- read.newick(url("https://raw.githubusercontent.com/Ensembl/ensembl-compara/release/97/scripts/pipeline/species_tree.vertebrates.branch_len.nw"))
  tree$edge.length <- NULL
  
  for (gene.name in gene.list) {
    print(gene.name)
    setwd(input.directory)
    summary.ok <- data.frame("gene.name" = character(), stringsAsFactors = FALSE)
    summary.error <-  data.frame("gene.name" = character(), stringsAsFactors = FALSE)
    
    #cambiar en vez de alignment type que diga suffix
    
    #  if (alignment.type == "omm_macse") {
    # The alignments drawn from omm_macse have the following structure:
    # suffix <- "_final_align_NT.aln"
    
    # The alignments drawn from prank/trimal have the following structure:
    # suffix <- "_nuc_trimmed.fas"
    
    setwd(input.directory)
    gene.name.omm <- paste0(gene.name, suffix)      
    dna <- readDNAStringSet(gene.name.omm)
    tmp <- attributes(dna)$ranges
    species <- names(tmp)
    species
    
    pruned.tree<-drop.tip(tree, tree$tip.label[-na.omit(match(species, tree$tip.label))])
    write.tree(pruned.tree)
    check <- length(species) == length(pruned.tree$tip.label)     #Sanity check species numb = tree tips numb = true
    check  
    unrooted.pruned.tree <- unroot(pruned.tree)
    destination.nwk <- paste(gene.name, '_pruned.tre', sep = "")
    file.destination.name <- paste0(tree.output.directory, destination.nwk)
    
    setwd(tree.output.directory)
    write.tree(unrooted.pruned.tree, file = destination.nwk)
    
    setwd(input.directory)
    
    width.seq <- dna@ranges@width[1]
    length.seq <- length(species)
    width.seq
    length.seq
    
    if (width.seq %% 3 == 0) { # width.seq %% 3 means triplets ok
      table.gene.name <- read.table(paste(gene.name.omm, sep = ""))
      max.index <- nrow(table.gene.name)+1
      temp <- data.frame(matrix(ncol=1, nrow=max.index))
      temp[1,1] <- as.character(paste("", length.seq, width.seq, sep = "\t"))
      temp[c(2:max.index),1] <- as.character(table.gene.name[,1])
      output.name <- paste(gene.name, ".phy", sep = "")
      phylip.file <- paste0(tree.output.directory,output.name)
      write.table(temp, phylip.file, col.names = FALSE, row.names = FALSE, quote = FALSE, eol = "\n")
      
      print ("OK")
      print(output.name)
      
      summary.gene.ok <-  data.frame("gene.name" = as.character(gene.name),stringsAsFactors = FALSE)
      summary.ok <- rbind(summary.ok,summary.gene.ok) 
      
      # human is our reference so it should be atop: open phylip.file and rearrange human at top

      
    } else {
      print ("Error: Total alignment length not %% 3")
      print(output.name)
      summary.gene.error <-  data.frame("gene.name" = as.character(gene.name),stringsAsFactors = FALSE)
      summary.error <- rbind(summary.error,summary.gene.error) 
    }
  }
 
      
      
     
     datestamp <- date()
  write.table(summary.ok, file = paste0(datestamp," Genes_OK.txt"), row.names = FALSE)
  write.table(summary.error, file = paste0(datestamp," Genes_Error.txt"), row.names = FALSE)
} 