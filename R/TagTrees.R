##Tag trees in foreground lineages

extract.tag.tree <- function(gene.list, m0.directory, classifier.directory, Ha.directory, Ho.directory) 
{
  setwd(m0.directory)
  for (gene.name in gene.list) {
    #extract tree
    #gene.name <- "RBFOX2"
    lines <- readLines(paste0("mlc_M0_", gene.name))
    lines.start.array <- startsWith(lines, "tree length =")
    pos.start <- which(lines.start.array == TRUE)
    lines.end.array <- startsWith(lines, "Detailed output identifying parameters")
    pos.end <- which(lines.end.array == TRUE)
    lineas.utiles <- lines[pos.start:pos.end]
    pos.enters <- which(lineas.utiles == "")
    utiles2 <- lineas.utiles[pos.enters[2]:pos.enters[3]]
    utiles2 <- paste(utiles2, sep = "", collapse = "")
    utiles2 <- trimws(utiles2)
    #tag tree in foreground
    
    tr <- read.tree(text = utiles2)
    plot(tr, show.node.label = TRUE)
    tr3 <- tr
    
    all.spp <- tr$tip.label
    
    setwd(classifier.directory)
    output.classified.gene <- paste0(gene.name, "_classification.csv")
    csv <- read.csv(output.classified.gene, sep = ";")
    colnames(csv) <- c("names", "class")
    csv[,"names"] <- tolower(csv[,"names"])
    csv <- csv[,c("names", "class")]
    
    fore.spp <- c()
    
    for  (species in all.spp) {
      
      
      if (csv[which(csv$names == species), 2] == "foreground") {
        fore.spp <-  c(fore.spp, species)      
      }
      
    }
    
    if (length(fore.spp) == 1) {
     
      pattern <- paste0(fore.spp,": ", "[:digit:]", ".", "[:digit:]+")
     # pattern <- paste0("homo_sapiens: ", "[:digit:]", ".", "[:digit:]+")
      str_view(utiles2, pattern)
      to.tag <- str_extract(utiles2, pattern)
      replacement <- paste0(to.tag, " #1")
      utiles3 <- str_replace(utiles2, to.tag, replacement)
      setwd(m0.directory)
      output.name <- paste0("tagged_mlc_M0_", gene.name, ".txt")
      writeLines(utiles3, con = output.name)
      file.copy(output.name, Ha.directory)
      file.copy(output.name, Ho.directory)
      file.copy(paste0(gene.name, ".phy"), Ha.directory)
      file.copy(paste0(gene.name, ".phy"), Ho.directory)
    } else {
      tr3 <- makeNodeLabel(tr3, "u", nodeList = list("#1" = fore.spp))
      plot(tr3, show.node.label = TRUE)
      setwd(m0.directory)
      output.name <- paste0("tagged_mlc_M0_", gene.name, ".txt")
      write.tree(tr3, file = output.name)
      
      file.copy(output.name, Ha.directory)
      file.copy(output.name, Ho.directory)
      file.copy(paste0(gene.name, ".phy"), Ha.directory)
      file.copy(paste0(gene.name, ".phy"), Ho.directory)
    }
    
    
   
    
   # file.remove(output.name)
  #  file.remove(paste0(gene.name, ".phy"))
    
    
  }
  
  
}