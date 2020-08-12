
##RUN M0 model
run.m0.paml <- function (gene.list, input.directory, output.directory, software.directory)
{
  
  
  #https://rest.ensembl.org/documentation/info/taxonomy_classification  tengo q mejorarlo con esto
  
  #
  codeml <- paste0(software.directory, "/codeml")
  ctl.original <- paste0(software.directory,"/M0.ctl.original")
  file.copy(codeml, output.directory)
  file.copy(ctl.original, output.directory)
  
  
  #gene.name <- "RBFOX1"
  for (gene.name in gene.list) {
    
    setwd(input.directory)
    
    nwk.file <- paste0(input.directory, gene.name, '_pruned.tre')
    gene.phylip <- paste0(input.directory, gene.name, ".phy")
    
    #set up everything locally for paml to work
    
    file.copy(nwk.file, output.directory)
    file.copy(gene.phylip,output.directory)
    
    setwd(output.directory)
    
    #m0 customization
    control.file.m0 <- readLines("m0.ctl.original")
    custom.control.file.m0 <- str_replace(control.file.m0,"pepito",gene.name)
    
    custom.control.file.m0.name <- paste0(gene.name,"_m0_ctrlfile")
    
    writeLines(custom.control.file.m0, custom.control.file.m0.name)
    
    
    #run m0 model
    
    setwd(output.directory)
    codeml.m0 <- paste0("./codeml ", custom.control.file.m0.name)
    system(codeml.m0)
    
    
  }
  
}
