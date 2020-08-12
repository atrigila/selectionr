#!/usr/bin/env Rscript
#Slow foreground classifier 
species.foreground.classifier <- function(gene.list, alignment.type, Taxon.id = "Mammalia", input.directory, output.directory)
{
  
  
  
  library(httr)
  library(jsonlite)
  library(xml2)
  library(stringr)
  #gene.list <- "MET"  
  
  setwd(input.directory)
  
  for (gene.name in gene.list) {
    
    #gene.name <- "MET"
    #Taxon.id <- "Mammalia"
    
    #gene.name <- c("procesadoRBFOX1")
    
    setwd(input.directory)
    
    if (alignment.type == "prank") {
      gen.start <- gene.name
      gen.end <- "_pranked.best.nuc_hmm.fasta"
      gen.file <- paste0(gene.name, gen.end)
      dna <- readDNAStringSet(gen.file)
      tmp <- attributes(dna)$ranges
      species <- names(tmp)
      species
    }
    if (alignment.type == "omm_macse") {
      gen.start <- gene.name
      gen.end <- "_final_align_NT.aln"
      gen.file <- paste0(gene.name, gen.end)
      dna <- readDNAStringSet(gen.file)
      tmp <- attributes(dna)$ranges
      species <- names(tmp)
      species
      print(gene.name)
    }
    
    
    classification.spp <- data.frame("names" = character(), "class" = character())
    
    for (specie in species) 
    {
      
      #specie<- "homo_sapiens"
      specie.1 <- word(string = specie, start = 1, end = 1, sep = "_")
      server <- "https://rest.ensembl.org/taxonomy/classification/"
      individual.query <- paste0(server,specie.1, "?")
      print(individual.query)
      
      Sys.sleep(1)
      r <- GET(paste(individual.query, sep = ""), content_type("application/json"))
      
      stop_for_status(r)
      
      # use this if you get a simple nested list back, otherwise inspect its structure
      # head(data.frame(t(sapply(content(r),c))))
      pepe <- (fromJSON(toJSON(content(r))))
      
      
      #Taxon.id <- "Mammalia"
      if(Taxon.id %in% pepe$scientific_name == "TRUE"){
        target.taxon.id.species <- data.frame("names" = specie, "class" = "foreground")
        classification.spp <- rbind(target.taxon.id.species, classification.spp)
      } else {
        target.taxon.id.species <- data.frame("names" = specie, "class" = "background")
        classification.spp <- rbind(target.taxon.id.species, classification.spp)
        
      }
      output.name <- paste0(gene.name, "_classification.csv")
      setwd(output.directory)
      write.csv2(classification.spp, output.name, row.names = FALSE)
      
    }
    
    
  }
  
  
}