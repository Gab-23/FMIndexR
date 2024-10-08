% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FM_index_from_FASTA.R
\name{FM_index_from_FASTA}
\alias{FM_index_from_FASTA}
\title{FM_index_from_FASTA}
\usage{
FM_index_from_FASTA(input, output, save = TRUE, compress = FALSE)
}
\arguments{
\item{input}{path to the FASTA file}

\item{output}{path to the directory in which to save the data structures}

\item{save}{default = TRUE, can be modified
to avoid saving data in the output directory}

\item{compress}{default = FALSE, can be modified
to compress the FM index, by downsampling the suffix array}
}
\value{
List-like object of type FM_index
}
\description{
The function generates an FM_index structure starting
from an input FASTA file.
The FM Index contains the suffix array (SA),
the Burrows-Wheeler Transform (BWT),
the Occurrencies matrix (Occ) and the Count array (C).
}
\details{
The FASTA file is parsed using a Biostrings function (readDNAStringSet)
in order to ensure a robust and reliable parsing method.

The FM index is compressed by downsampling the suffix array,
taking one every 32 suffixes.
}
\examples{
# Example creation of an FM index from a FASTA file
input_file <- system.file("extdata", "prova.txt", package = "FMIndexR")
output_path <- system.file("output", package = "FMIndexR")

FM_index <- FM_index_from_FASTA(input_file,
                                output_path,
                                save = FALSE,
                                compress = FALSE)
FM_index
}
