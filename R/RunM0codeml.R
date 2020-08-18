#' Run PAML M0 model
#'
#' @param gene.list Gene name or file name
#' @param input.directory Directory where newick and phylip files are located
#' @param output.directory Directory where M0 optimized files will be written
#' @param software.directory Directory where PAML software is locateed
#'
#' @importFrom stringr str_replace
#'
#' @return A mlc file along with txt others files written by PAML. We mostly
#'  make use of the mlc file with the optimized tree.
#' @export run.m0.paml

run.m0.paml <- function(gene.list, input.directory, output.directory, software.directory) {

  codeml <- paste0(software.directory, "/codeml")
  ctl.original <- paste0(software.directory,"/M0.ctl.original")
  file.copy(codeml, output.directory)
  file.copy(ctl.original, output.directory)


  for (gene.name in gene.list) {

    setwd(input.directory)

    nwk.file <- paste0(input.directory, gene.name, '_pruned.tre')
    gene.phylip <- paste0(input.directory, gene.name, ".phy")

    # Set up everything locally for paml to work

    file.copy(nwk.file, output.directory)
    file.copy(gene.phylip,output.directory)

    setwd(output.directory)

    # M0 customization
    control.file.m0 <- readLines("m0.ctl.original")
    custom.control.file.m0 <- stringr::str_replace(control.file.m0,"pepito",gene.name)

    custom.control.file.m0.name <- paste0(gene.name,"_m0_ctrlfile")

    writeLines(custom.control.file.m0, custom.control.file.m0.name)


    # Run M0 model

    setwd(output.directory)
    codeml.m0 <- paste0("./codeml ", custom.control.file.m0.name)
    system(codeml.m0)


  }

}
