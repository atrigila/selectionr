### PRANK Aligner

# Align sequences with prank
#gene.name <- gene.list[1]

prank.align <- function (gene.list, input.file.directory, output.file.directory, software.directory) {
  setwd(output.file.directory) 
  file.copy(paste0(software.directory, "/prank"), output.file.directory)
  # Another option: use species tree phylogeny
  # file.copy(paste0(software.directory, "EnsemblTree.txt"), deep.prank.directory) #Phylogenetic guide tree for prank
  # more complicated:   
  # prank.align <- paste0("./prank ","-d=", "'", custom.gene,"'"," -t=EnsemblTree.txt", " -o=" , gene.name , "_pranked" , " -translate", " -F", " -prunetree -once")
  
  file.copy(paste0(software.directory, "/faTrans"), output.file.directory)
  
  
  
  for (gene.name in gene.list) {
    
    file.type <- '.fasta'
    custom.gene <- paste(input.file.directory, gene.name, file.type, sep = "")
    print(custom.gene)
    

    if (file.exists(custom.gene) == TRUE) {
      setwd(output.file.directory) 
      
      #command.translate <- paste0("./faTrans '", custom.gene,"' " ,gene.name, "_AA.fasta")
      #system(command.translate)
      
      #prank.align <- paste0("./prank ","-d=", "'", gene.name,"_AA.fasta'", " -o=" , gene.name , "_pranked") # simplest alignment in aminoacids
      #system(prank.align)
      
      prank.align <- paste0("./prank ","-d=", "'", custom.gene,"'"," -o=" , gene.name , "_pranked" , " -translate", " -F")
      system(prank.align)
      
      
      
    } else {
      print("Error alignining gene")
    }
  }
  setwd(output.file.directory) 
  file.remove(paste0(output.file.directory, "/prank"))
  file.remove(paste0(output.file.directory, "/faTrans"))
  
#  file.remove(paste0(software.directory, "EnsemblTree.txt"))
  
  
}
