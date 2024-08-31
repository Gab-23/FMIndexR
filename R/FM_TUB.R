#' FM_TUB
#'
#' An example FM index created from a FASTA file.
#'
#' The FASTA file is stored in inst/extdata/TUB_FASTA.txt
#'
#' @docType data
#' @usage data(FM_TUB)
#' @format An object of class FM_index, containing:
#' \describe{
#'   \item{Suffix Array}{The Suffix Array of the sequence}
#'   \item{BWT}{The Burrows-Wheeler Transform (BWT) of the sequence}
#'   \item{Occ matrix}{The Occurrencies matrix of the BWT}
#'   \item{Count array}{The Count array of the sequence}
#'   \item{FASTA header}{The FASTA header of the sequence}
#' }
#' @examples
#' data(FM_TUB)
#' @source \url{https://www.ncbi.nlm.nih.gov/nuccore/NM_177972.3}
"FM_TUB"
