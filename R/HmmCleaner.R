## HMM Clean in AA, TransferClean to NT and output

HmmClean <- function(gene.list, suffix_AA = "_pranked.best.pep.fas", suffix_NT = "_pranked.best.nuc.fas", input.align.directory, output.clean.directory) {
   
   library(tidyverse)
  
  summary.transfer.clean <- data.frame("Gene.name" = character()) 
  
    for (gene.name in gene.list) {
       
   #   gene.alignment <- paste0(gene.name,"_pranked.best.pep.fas")
      gene.alignment <- paste0(gene.name, suffix_AA)
      hmm.clean <- paste0("HmmCleaner.pl '", gene.alignment,"'")
      setwd(input.align.directory)
      system(hmm.clean) # https://metacpan.org/pod/distribution/Bio-MUST-Apps-HmmCleaner/bin/HmmCleaner.pl
      
      
      gene.nt.align <- paste0(gene.name, suffix_NT)
      gene.aa.log <- paste0(gene.name, str_remove(suffix_AA, ".fas"), "_hmm.log")
      transfer.clean <- paste0("transferCleaner.pl '", gene.nt.align,"'", " -log='", gene.aa.log,"'")
      setwd(input.align.directory)
      system(transfer.clean) # https://metacpan.org/pod/distribution/Bio-MUST-Apps-HmmCleaner/bin/transferCleaner.pl
      
      copy.clean.aa <-  list.files(path = input.align.directory, pattern = "_hmm")
      file.copy(copy.clean.aa, output.clean.directory)
      file.remove(copy.clean.aa)
      
      copy.transfer.nt <- list.files(path = input.align.directory, pattern = "_cleaned.ali")
      file.copy(copy.transfer.nt, output.clean.directory)
      file.remove(copy.transfer.nt)
      
      setwd(output.clean.directory)
      
      name.cleaned <- paste0(gene.name, str_remove(suffix_NT, ".fas"),"_cleaned.ali")
      if (file.exists(name.cleaned) == TRUE) {
        cleaned.ali <- readLines(name.cleaned) 
        cleaned.ali <- str_remove(cleaned.ali, "#")
        cleaned.ali <- str_replace_all(cleaned.ali, pattern = "\\*", replacement =  "\\-")
        cleaned.ali <- str_replace(cleaned.ali, pattern = "\\t\\t>", replacement = "")
        
        write.table( x =cleaned.ali , file = paste0(gene.name, "_cleaned.fasta"), quote = FALSE, row.names = FALSE, col.names = FALSE)
        
      } else {
        print("Not transferCleaned")
        gene.id <- data.frame("Gene.name" = gene.name)
        summary.transfer.clean <- rbind(gene.id, summary.transfer.clean)
        setwd(output.clean.directory)
        file.remove(copy.transfer.nt)
        file.remove(copy.clean.aa)
        
        
      }

    
    }
  write.table(summary.transfer.clean, file = "1.NotTransferCleaned.txt")
    
}


