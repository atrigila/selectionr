### Foreground classifier 

# Foreground classifier (medium version):

# Extract all unique species present in all alignments

species.classifier <- function(gene.list, suffix ="_nuc_trimmed.fas" , input.directory, Taxon.id = "Mammalia", classification.output.directory) {
  library(jsonlite)

all.species <- c()
  for (gene.name in gene.list) {
    gene.dna <- paste0(input.directory, gene.name,suffix)
    dna <- readDNAStringSet(gene.dna)
    species.per.gene <- names(dna)
    all.species <- c(all.species,species.per.gene)
  }

  species <- unique(all.species)
  
  
# Classify all unique species if they are foreground or background 
# Ask Ensembl if the species genus belongs to my target foreground taxon id
  classification.spp <- data.frame("names" = character(), "class" = character())
  
  for (specie in species)  {
    
    #specie<- "homo_sapiens"
    specie.1 <- word(string = specie, start = 1, end = 1, sep = "_")
    server <- "https://rest.ensembl.org/taxonomy/classification/"
    individual.query <- paste0(server,specie.1, "?")
    print(individual.query)
    
    Sys.sleep(1)
    r <- GET(paste(individual.query, sep = ""), content_type("application/json"))
    
    stop_for_status(r)
    Sys.sleep(1)
    pepe <- (fromJSON(toJSON(content(r))))
    
    
    #Taxon.id <- "Mammalia"
    if(Taxon.id %in% pepe$scientific_name == "TRUE"){
      target.taxon.id.species <- data.frame("names" = specie, "class" = "foreground")
      classification.spp <- rbind(target.taxon.id.species, classification.spp)
    } else {
      target.taxon.id.species <- data.frame("names" = specie, "class" = "background")
      classification.spp <- rbind(target.taxon.id.species, classification.spp)
      
    }
    
     
  }
  

for (gene.name in gene.list) {
  
  gene.dna <- paste0(input.directory, gene.name,suffix)
  dna <- readDNAStringSet(gene.dna)
  species.per.gene <- names(dna)
  
  classification.spp.out <- data.frame("names" = character(), "class" = character())
  
  foreground.species <- classification.spp[classification.spp$class == "foreground", "names"]
  

  for (specie in species) {
    if(specie %in% foreground.species == "TRUE"){
      target.taxon.id.species <- data.frame("names" = specie, "class" = "foreground")
      classification.spp.out <- rbind(target.taxon.id.species, classification.spp.out)
    } else {
      target.taxon.id.species <- data.frame("names" = specie, "class" = "background")
      classification.spp.out <- rbind(target.taxon.id.species, classification.spp.out)
      
    } 
    output.name <- paste0(classification.output.directory, gene.name, "_classification.csv")
    write.csv2(classification.spp.out, output.name, row.names = FALSE)
    
  }
  
  
}
  
}