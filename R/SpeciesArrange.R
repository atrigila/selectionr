### Function to make a species be a top of the alignment
# Useful when some type of aligners rearrange species order
# Works for dna

species.top <-  function(gene.list, suffix = "_nuc_trimmed.fas", species.atop = "homo_sapiens", input.directory, output.directory) {
  library(Biostrings)
  for (gene.name in gene.list) {
    #example suffix <- "_nuc_trimmed.fas"
    gene.file <- paste0(input.directory, gene.name, suffix)
    dna <- readDNAStringSet(gene.file)
    
    moveToTop <- function(x, n) c(x[n], x[-n]) # https://support.bioconductor.org/p/127125/
    
    # which species atop? species.atop <- "homo_sapiens"
    index <- which(species.atop == names(dna))
    new.order <- moveToTop(dna, index)
    
    output.name <- paste0(gene.name,"_ordered" ,suffix)
    
    writeXStringSet(x= new.order, filepath =  paste0(output.directory, output.name))
    
  }
  
  
}