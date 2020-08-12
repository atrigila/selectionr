

### Get gene original sequence, without gaps and without masked nucleotides

#input.directory <- "/Users/Usuario/Desktop/selection/results extensive/results omm_macse_extensive/2. Align omm_macse_extensive/"
#gene.name <- 'ADGRB1'
#suffix <- '_final_align_NT.aln'
#output.directory <- "/Users/Usuario/Desktop/"
#maskFullFasta <- "/Users/Usuario/Desktop/selection/results extensive/results omm_macse_extensive/2. Align omm_macse_extensive/ADGRB1/ADGRB1_maskFull_detail.fasta"
#beb.sites.directory <- "/Users/Usuario/Desktop/selection/results extensive/results omm_macse_extensive/9.Summary omm_macse_extensive/"
positions.provean <- function(gene.name,
                              input.directory = "/Users/Usuario/Desktop/selection/results extensive/results omm_macse_extensive/2. Align omm_macse_extensive/" ,
                              suffix = '_final_align_NT.aln',
                              output.directory ,
                              classifier.directory = "/Users/Usuario/Desktop/selection/results extensive/results omm_macse_extensive/5. Classify omm_macse_extensive/",
                              beb.sites.directory = "/Users/Usuario/Desktop/selection/results extensive/results omm_macse_extensive/9.Summary omm_macse_extensive/") {



  # 0.1.  Set libraries and useful dataframes
  library(stringr)
  library(seqinr)
  library(Biostrings)
  library(jsonlite)

  print(gene.name)



  beb.sites.human.newnumber <- data.frame("Gene Name" = integer(),
                                          "ENSEMBL Protein ID" = character(),
                                          "Position number" = integer(),
                                          "Human BEBs Letter" = character(),
                                       #   "Outgroup BEB name" =character(),
                                          "Outgroup BEBs" = character(),
                                          "New Position number (nogaps)" = integer(),
                                          "New position BEBs letter check" = character(),
                                          "Original number" = integer(),
                                          "Original letter check" = character())


  # 0.2. Get name of ENSEMBL ID to which original letter corresponds to

  # This is slow but and there is probably a more efficient way to query the database. #

  target_taxon <-9443 # This should be adapted to the set of species in use. It is ok in my case for mammals as I only need info for human.
  url.head <- 'https://jul2019.rest.ensembl.org/homology/symbol/human/'
  url.taxon <-paste('?target_taxon=',target_taxon,';', sep="")
  url.tail <- 'content-type=application/json;type=orthologues'
  destfile.custom <- paste(gene.name, '.json', sep="")
  url.custom <- paste(url.head, gene.name, url.taxon, url.tail, sep="")
  print(url.custom)

  JSON.object <- tryCatch(fromJSON(url.custom),
                          error = function (e) {
                            print(paste("This gene does not exist:", gene.name))
                            return(NULL)})

  ensembl.protein.name<- JSON.object$data$homologies[[1]]$source$protein_id[1]

  Sys.sleep(1)

  # 1. Load alignment used for PAML codeml.  Remove gaps.

  final.alignment <- paste0(input.directory,"/",gene.name,"/", gene.name, suffix)

  final.align <- readDNAStringSet(final.alignment)
  number.human <- which(names(final.align) == "homo_sapiens")
  second.name <- names(final.align[number.human + 1])

  final.align <- readLines(final.alignment)
  human.final.align.start <- startsWith(final.align, ">homo_sapiens") #can be modified to accept any name
  pos.start.final.align <- which(human.final.align.start == TRUE)

  human.final.align.end <- endsWith(final.align, paste0(">",second.name))
  pos.end.final.align <- which(human.final.align.end == TRUE)

  useful.lines.final.align <- str_remove(final.align[pos.start.final.align:pos.end.final.align], paste0(">",second.name))

  split.useful.lines <- strsplit(useful.lines.final.align[2],"")
  letters <- split.useful.lines[[1]]
  index.gaps <- which(str_detect(split.useful.lines[[1]], "-"))

  sequence.minus.gaps <- letters[-index.gaps]

  length(sequence.minus.gaps)

  # 2.  Read the original masked sequence, select only the human sequence (tracing.file)

  maskFullFasta <- paste0(input.directory, "/", gene.name, "/", gene.name, "_maskFull_detail.fasta")
  s <- readDNAStringSet(maskFullFasta)
  second.name <- names(s)[2]

  maskFulldetail <- readLines(maskFullFasta)
  human.maskfulldetail.start <- startsWith(maskFulldetail, ">homo_sapiens")
  pos.start <- which(human.maskfulldetail.start == TRUE)

  human.maskfulldetail.end <- endsWith(maskFulldetail, paste0(">",second.name))
  pos.end <- which(human.maskfulldetail.end == TRUE)

  useful.lines <- str_remove(maskFulldetail[pos.start:pos.end], paste0(">",second.name))

  split.useful.lines <- strsplit(useful.lines[2],"")
  letters.masked <- split.useful.lines[[1]]
  index.lowecase <- which(str_detect(letters.masked, "^[:lower:]+$"))

  sequence.minus.masked.nucleotides <- letters.masked[-index.lowecase]

  length(sequence.minus.masked.nucleotides)

  print("Are ungapped and unmasked sequences the same?")
  length(sequence.minus.gaps) == length(sequence.minus.masked.nucleotides) #sanity check

  # 3. For both sequences (ungapped and unmasked), get their correspoding protein translation

  translate.gapped <- seqinr::translate(letters)
  translate.gapped <- str_replace_all(translate.gapped, pattern = "X", replacement = "-")

  translate.ungapped <- seqinr::translate(sequence.minus.gaps)
  translate.unmasked <- seqinr::translate(sequence.minus.masked.nucleotides)

  print("Are there any differences between ungapped and unmasked sequences?")
  #isFALSE(translate.ungapped == translate.unmasked) # sanity check, should be false!

  # 4. Read BEB sites and find them in the alignment used for paml (letters)

  beb.sites.file <- read.delim(paste0(beb.sites.directory, "BEBSites.txt"))
  sites.for.gene <- beb.sites.file[which(beb.sites.file$gene.name == gene.name),"BEBsites"]

  pattern <- "[:digit:]+[:space:]"
  sites.ok <- str_extract_all(sites.for.gene, pattern)
  sites.ok <- as.integer(unlist(sites.ok))
  sites.ok <- str_replace(sites.ok, " ", ",")
  sites.ok <- as.integer(sites.ok)

  # Get the outgroup seq

  classification <- paste0(classifier.directory, gene.name, "_classification.csv")
  classification <- read.csv(classification, sep= ';', stringsAsFactors = FALSE)

  index.backgrounds <- which(classification$class == 'background')
  background.names <- classification$names[index.backgrounds]

  msa.seq <-seqinr::read.fasta(final.alignment)
  background.names <- intersect(names(msa.seq),background.names)

  msa.seq.background <- msa.seq[background.names]#get aminoacids for a background outgroup sequences
  outgroup.seq.aa <- lapply(msa.seq.background, seqinr::translate)


  # 6. Get the original position in the no-gaps, no-mask sequence: letters.masked

  # Translate the letters.masked taking care of the masked and unmasked sequences

  replaced.lowercase <- str_replace_all(letters.masked,pattern = "^[:lower:]+$", "-")
  original.replaced.masked <- seqinr::translate(replaced.lowercase)
  original.replaced.masked <- str_replace_all(original.replaced.masked, pattern = "X", replacement = "-")



  # 5. Read BEB sites and find them in the original alignment (translate.gapped)
  #number <- 1371
  for (number in sites.ok){
    how.many.gaps.to.site <- sum(str_count(translate.gapped[1:number], "-")) # Count how many gaps until the site number.
    new.number <- number-how.many.gaps.to.site # This should be equal to the tracing file
    translate.ungapped[new.number] == translate.unmasked[new.number] # sanity check

    gaps <- which(original.replaced.masked == "-")
    inicio <- c(gaps[1]) #el inicio siempre es la primera posicion donde hay un gap
    fin <- c()

    indices <- 1:(length(gaps)-1)

    for (i in indices) {
      if (gaps[i+1] != (gaps[i])+1){
        inicio <- c(inicio, gaps[i+1])
        fin <- c(fin, gaps[i])
      }
    }

    fin <- c(fin,gaps[length(gaps)])
    inicios.y.fines <- data.frame(inicio,fin)
    inicios.y.fines["length"]<- NA
    inicios.y.fines$length <- inicios.y.fines$fin - inicios.y.fines$inicio + rep(1,nrow(inicios.y.fines))

    potenciales.posiciones <- which(translate.ungapped[new.number] == original.replaced.masked)

    #potencial.posicion <- 180
    #bloque.gap <- 1
 ###### TEST ########




    for (potencial.posicion in potenciales.posiciones) {
      how.many.gaps.to.potencial.posicion <- sum(str_count(original.replaced.masked[1:potencial.posicion], "-"))
      if(potencial.posicion-how.many.gaps.to.potencial.posicion == new.number){
        original.number <- potencial.posicion
        original.letter <- original.replaced.masked[potencial.posicion]
      }
    }


  #  number <- 107
    aminoacids <- c()
    for (name in background.names) {
      letter.outgroup <- outgroup.seq.aa[[name]][[number]]
      aminoacids <- c(aminoacids,letter.outgroup)
    }
    aminoacids <- aminoacids[!aminoacids %in% "X"]


    calculate_mode <- function(x) {
      uniqx <- unique(na.omit(x))
      uniqx[which.max(tabulate(match(x, uniqx)))]}


    most.frequent.outgroup.aminoacid <- calculate_mode(aminoacids)



######################
    # 7. Output results



    single.site.gene <- data.frame("Gene Name" = as.character(gene.name),
                                   "ENSEMBL Protein ID" = as.character(ensembl.protein.name),
                                   "Position number" = number,
                                   "Human BEBs Letter" = translate.gapped[number],
                         #          "Outgroup BEB name" =as.character(first.background),
                                   "Outgroup BEBs" = most.frequent.outgroup.aminoacid,
                                   "New Position number (nogaps)" = new.number,
                                   "New position BEBs letter check" = translate.unmasked[new.number],
                                   "Original number" = original.number,
                                   "Original letter check" = original.letter,
                                    stringsAsFactors = FALSE)

    beb.sites.human.newnumber <- rbind(beb.sites.human.newnumber,single.site.gene)

  }



  # 8. Write to file

  to.provean.sites <- paste(beb.sites.human.newnumber$Original.number,
                            beb.sites.human.newnumber$Original.letter.check,
                            beb.sites.human.newnumber$Outgroup.BEBs, sep=",")

  fileConn<-file(paste0(output.directory, gene.name, "_PROVEAN_OriginalPositions.txt"))
  writeLines(to.provean.sites, fileConn)
  close(fileConn)


  name.human<- paste0(output.directory, gene.name, "_human_AA.fas")
  write.fasta(seqinr::translate(letters.masked), names = "homo_sapiens", file.out =  name.human)

  output.table <- beb.sites.human.newnumber[,c("Gene.Name", "ENSEMBL.Protein.ID","Original.number","Original.letter.check","Outgroup.BEBs")]
  write.table(output.table, file = paste0(output.directory, gene.name, "_info_id.txt"), sep= "\t", row.names= F, quote = F)
  }





