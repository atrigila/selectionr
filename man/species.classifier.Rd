% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SpeciesClassifier.R
\name{species.classifier}
\alias{species.classifier}
\title{Foreground classifier}
\usage{
species.classifier(
  gene.list,
  suffix = "_nuc_trimmed.fas",
  input.directory,
  Taxon.id = "Mammalia",
  classification.output.directory
)
}
\arguments{
\item{gene.list}{character; Gene list (file name) of MSAs to check for species names}

\item{suffix}{character; Suffix of the files (such as: ".fasta")}

\item{input.directory}{character; Directory where the files for the MSA are located}

\item{Taxon.id}{character; Foreground taxon to check if the species belongs to}

\item{classification.output.directory}{character; Directory where the csv will be written}
}
\value{
A csv with species names and class (either foreground or background)
}
\description{
Extract all unique species present in all alignments and assign to
foreground and background sets if they belong to a certain taxon id or not
}
