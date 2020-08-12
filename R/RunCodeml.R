run.paml <- function (gene.list, software.directory, Ho.directory, Ha.directory)
{
  
  file.copy(paste0(software.directory, "/branch-site_Ho.ctl"), Ho.directory)
  file.copy(paste0(software.directory, "/branch-site_Ha.ctl"), Ha.directory)
  file.copy(paste0(software.directory, "/codeml"), Ho.directory)
  file.copy(paste0(software.directory, "/codeml"), Ha.directory)
  
  
  library(tidyverse)
  library(stringr)
  
  #gene.name<-"COL11A1"
  for (gene.name in gene.list)
  {
    gene.phylip <- paste(gene.name, ".phy", sep = "")
    
    #Ho customization
    setwd(Ho.directory)
    control.file.ho <- readLines("branch-site_Ho.ctl")
    custom.control.file.ho <- str_replace(control.file.ho,"pepito",gene.name)
    
    custom.control.file.ho.name <- paste0(gene.name, "_branch-site_Ho.ctl")
    print(custom.control.file.ho.name)
    
    writeLines(custom.control.file.ho, custom.control.file.ho.name)
    
    codeml.ho <- paste0("./codeml ", custom.control.file.ho.name)
    system(codeml.ho)
    
    
    #Ha customization
    setwd(Ha.directory)
    control.file.ha <- readLines("branch-site_Ha.ctl")
    custom.control.file.ha <- str_replace(control.file.ha,"pepito",gene.name)
    custom.control.file.ha.name <- paste0(gene.name, "_branch-site_Ha.ctl")
    print(custom.control.file.ha.name)
    
    
    writeLines(custom.control.file.ha, custom.control.file.ha.name)
    
    
    codeml.ha <- paste0("./codeml ", custom.control.file.ha.name)
    
    system(codeml.ha)
    
  }
}
