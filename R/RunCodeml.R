#' Run Paml Codeml
#'
#' @param gene.list A gene list or file name
#' @param software.directory Directory where PAML Codeml executables are
#' @param Ho.directory Directory where the null hypothesis files will be written
#' @param Ha.directory Directory where the alternative hypothesis will be written
#'
#' @importFrom stringr str_replace
#'
#' @return Several txt files after the codeml run (such as the mlc)
#' @export run.paml
#'
#'

run.paml <- function (gene.list, software.directory, Ho.directory, Ha.directory) {

  file.copy(paste0(software.directory, "/branch-site_Ho.ctl"), Ho.directory)
  file.copy(paste0(software.directory, "/branch-site_Ha.ctl"), Ha.directory)
  file.copy(paste0(software.directory, "/codeml"), Ho.directory)
  file.copy(paste0(software.directory, "/codeml"), Ha.directory)

  for (gene.name in gene.list)  {
    gene.phylip <- paste(gene.name, ".phy", sep = "")

    #Ho customization
    setwd(Ho.directory)
    control.file.ho <- readLines("branch-site_Ho.ctl")
    custom.control.file.ho <- stringr::str_replace(control.file.ho,"pepito",gene.name)

    custom.control.file.ho.name <- paste0(gene.name, "_branch-site_Ho.ctl")
    print(custom.control.file.ho.name)

    writeLines(custom.control.file.ho, custom.control.file.ho.name)

    codeml.ho <- paste0("./codeml ", custom.control.file.ho.name)
    system(codeml.ho)


    #Ha customization
    setwd(Ha.directory)
    control.file.ha <- readLines("branch-site_Ha.ctl")
    custom.control.file.ha <- stringr::str_replace(control.file.ha,"pepito",gene.name)
    custom.control.file.ha.name <- paste0(gene.name, "_branch-site_Ha.ctl")
    print(custom.control.file.ha.name)

    writeLines(custom.control.file.ha, custom.control.file.ha.name)

    codeml.ha <- paste0("./codeml ", custom.control.file.ha.name)
    system(codeml.ha)

  }
}
