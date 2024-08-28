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
#'   \item{Suffix Array}{Description of the Suffix Array}
#'   \item{BWT}{Description of the Burrows-Wheeler Transform (BWT)}
#'   \item{Occ matrix}{Description of the Occurrence matrix}
#'   \item{Count array}{Description of the Count array}
#'   \item{FASTA header}{Description of the FASTA header}
#' }
#' @examples
#' data(FM_TUB)
#' @source \url{https://www.ncbi.nlm.nih.gov/nuccore/NM_177972.3}
"FM_TUB"
