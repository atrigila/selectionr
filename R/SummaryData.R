#' Summary statistics and positively selected sites
#'
#' This function allows you to summarize several mlc files from PAML into a single table.
#'
#' @param gene.list Gene name (file name) of the MSA
#' @param Ho.directory Directory where PAML null hypothesis are located
#' @param Ha.directory Directory where PAML alternative hypothesis are located
#' @param summary.directory Directory summarized PAML files will be located
#' @keywords summary
#' @importFrom stringr str_extract
#' @import dplyr
#' @export summ.statistics.paml
#' @return A table with the selected sites for each gene

summ.statistics.paml <- function (gene.list, Ho.directory, Ha.directory, summary.directory) {


  summary.statistics.table <- data.frame(gene.name = character(), HaLnL = character(), HoLnL = character (), stringsAsFactors = FALSE)
  error.table <- data.frame(gene.name = character())

  for (gene.name in gene.list){
    print(gene.name) #gene.name <- "KITLG"
    # 1. Read H0 file
    file.name  <- paste0("mlc_H0_",gene.name)
    file.ho.location <- paste0(Ho.directory,"/", file.name)
    lines.ho <- readLines(file.ho.location)

    # 2. Read Ha file
    file.name  <- paste0("mlc_Ha_",gene.name)
    file.ha.location <- paste0(Ha.directory,"/", file.name)
    lines.ha <- readLines(file.ha.location)

    # 2.1. Which lines, start with "lnL(ntime:  0  np:  5): "?

    lines.start.array.ho <- startsWith(lines.ho, "lnL")
    pos.start.ho <- which(lines.start.array.ho == TRUE)

    lines.start.array.ha <- startsWith(lines.ha, "lnL")
    pos.start.ha <- which(lines.start.array.ha == TRUE)


    # 3. Which lines end with "+0.000000"?

    lines.end.array.ho <- endsWith(lines.ho, "+0.000000")
    pos.end.ho <- which(lines.end.array.ho == TRUE)

    lines.end.array.ha <- endsWith(lines.ha, "+0.000000")
    pos.end.ha <- which(lines.end.array.ha == TRUE)

    # 4. indices util inicio y final, lineas utiles <- lines[inicio, final]
    if(length(pos.start.ho & pos.start.ha) == 0) {
      print("ERROR")
      error.table.gene <- data.frame("gene.name" = as.character(gene.name))
      error.table <- rbind(error.table,error.table.gene)
    } else {


    lineas.utiles.ho <- lines.ho[pos.start.ho:pos.end.ho]
    clean.useful.lines.ho <- str_extract(lineas.utiles.ho,"-[:digit:]+\\.[:digit:]+")


    lineas.utiles.ha <- lines.ha[pos.start.ha:pos.end.ha]
    clean.useful.lines.ha <- str_extract(lineas.utiles.ha,"-[:digit:]+\\.[:digit:]+")
    options(digits = 11)
    clean.useful.lines.ha <- as.numeric(clean.useful.lines.ha)

    #5. Create a data frame with Ha Lnl, Ho Lnl a
    summary.statistics.gene <- data.frame("gene.name" = as.character(gene.name), "HaLnL" = clean.useful.lines.ha, "HoLnL" =  clean.useful.lines.ho, stringsAsFactors = FALSE)

    summary.statistics.table <- rbind(summary.statistics.table,summary.statistics.gene)
    } #fin del else
  }

  summary.statistics.table[,2] <- as.numeric(summary.statistics.table[,2])
  summary.statistics.table[,3] <- as.numeric(summary.statistics.table[,3])

  summary.output <- summary.statistics.table %>%
    mutate("LRT" = 2*(summary.statistics.table[,2]-summary.statistics.table[,3]))

  summary.output$p.value <- pchisq(summary.output$LRT, df=1, lower.tail = FALSE)


  setwd(summary.directory)
  write.table(summary.output, "Summary Statistics PAML.txt", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
  write.table(error.table, "Genes not processed - Error.txt", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

  ## Detect BEB sites

  table <- read.delim(paste0(summary.directory,"Summary Statistics PAML.txt"))
  table$BEBsites <- NA

  i <- 1
  for (i in 1:nrow(table)) {
    if (table$p.value[i] < 0.05) {
      # escanea el archivo y encontra la parte con los sitios beb, guarda aquellos q terminen en *


      file.name  <- paste0("mlc_Ha_",table[i,1])
      file.ha.location <- paste0(Ha.directory,"/", file.name)
      lines.ha <- readLines(file.ha.location)

      lines.start.array.ha <- startsWith(lines.ha, "Bayes Empirical Bayes")
      pos.start.ha <- which(lines.start.array.ha == TRUE)

      lines.end.array.ha <- endsWith(lines.ha, "The grid (see ternary graph for p0-p1)")
      pos.end.ha <- which(lines.end.array.ha == TRUE)

      lineas.utiles.ha <- lines.ha[pos.start.ha:pos.end.ha]
      clean.useful.lines.ha <- str_extract(lineas.utiles.ha, ("\\d+\\s[:upper:]\\s\\d+\\.\\d+\\*"))

      clean.useful.lines.ha<- clean.useful.lines.ha[!is.na(clean.useful.lines.ha)]
      clean.useful.lines.ha <- paste(clean.useful.lines.ha, collapse = ", ")
      if (clean.useful.lines.ha == "") {
        table$BEBsites[i] <- "No sites w>1"
      } else {
        table$BEBsites[i] <- clean.useful.lines.ha
      }
    } else {
      table$BEBsites[i] <- "NS"
    }
  }
  return(table)


getwd()
write.table(table, file = paste0(summary.directory,"BEBSites.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")



  }


