% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EnsemblDownload.R
\name{download.orthologues}
\alias{download.orthologues}
\title{Download orthologues from Ensembl}
\usage{
download.orthologues(
  server = "https://jan2020.rest.ensembl.org",
  gene.name,
  type = c("symbol", "Ensembl"),
  target_species = NULL,
  target_taxon
)
}
\arguments{
\item{server}{Version of Ensembl to use}

\item{gene.name}{character; gene name to download (as symbol or Ensembl gene ID)}

\item{type}{character; one of "gene symbol" or "Ensembl"}

\item{target_species}{character; Use with gene symbol only.}

\item{target_taxon}{numeric; Taxon to download species.}
}
\value{
character stringr; A fasta file with sequences
}
\description{
An improved version to download orthologues from Ensembl.
}
\examples{
\dontrun{
x <- c("ENSAMXG00005000412", "ENSAMXG00005000517")
ex.genes <- lapply(setNames(x,x), download.orthologues,
server = "https://jan2020.rest.ensembl.org", type = "Ensembl",
target_taxon = 7898)}


}
