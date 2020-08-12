### RUN TRIMAL

# Input aligned sequences as AA, trim them and then backtranslate to cds 

trimal.clean <- function (gene.list, suffix = "_pranked.best.fas", alignment.directory, trim.output.directory, unaligned.directory, software.directory) {
  
  # suffix will depend on the alignment method used
  # example: suffix <- "_pranked.best.fas"
  # gene.name <- gene.list[1]
  
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