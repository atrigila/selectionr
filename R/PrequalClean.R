#### Function to use PREQUAL cleaner

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