# Function to get genomic location of protein sites for ensembl proteins
#ensembl.protein.name <- "ENSP00000430945"
#original.number <- 1269

genomic.locator <- function(ensembl.protein.ID,original.number){
  library(httr)
  library(jsonlite)
  library(xml2)
  
  server <- "https://jul2019.rest.ensembl.org"
  ext <- paste0("/map/translation/", ensembl.protein.ID,"/")
  sites <- paste0(original.number,"..",original.number,"?")
  r <- GET(paste(server, ext, sites, sep = ""), content_type("application/json"))
  
  
  stop_for_status(r)
  
  Sys.sleep(1)
  
  str <- head(fromJSON(toJSON(content(r))))
  chr <- paste0("chr",str$mappings$seq_region_name)
  start <- str$mappings$start[[1]]
  end <- str$mappings$end[[1]]
  
  df <- data.frame("gene.name"  = (ensembl.protein.ID),
    "chr" = c(chr),
             "start" = c(start),
             "end" = c(end), stringsAsFactors=FALSE)
  
  return(df)

}

# Example usage
location <- "/Users/Usuario/Desktop/selection/results extensive/results omm_macse_extensive/10.Provean/"
all.genes.info <- read.delim("/Users/Usuario/Desktop/selection/results extensive/results omm_macse_extensive/10.Provean/All_Genes_info_id.txt")


#gene.name <- "ADGRB1"
locations <- data.frame("gene.name"  = character(), "chr" = character(), "start" = integer(), "end"= integer())

for (i in 1:nrow(all.genes.info)){
  gene.symbol <- all.genes.info$Gene.Name[i]
  gene.name <- all.genes.info$ENSEMBL.Protein.ID[i]
  original.number <- all.genes.info$Original.number[i]
  
  save.location <- genomic.locator(ensembl.protein.ID =  gene.name, original.number = original.number)
  locations <- rbind(locations, save.location)
}

all.info <- cbind.data.frame(all.genes.info,locations)

write.table(all.info, file = "/Users/Usuario/Desktop/selection/results extensive/results omm_macse_extensive/10.Provean/Genomic_locations.txt", sep = "\t", quote = F)


