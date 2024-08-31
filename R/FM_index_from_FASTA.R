#' FM_index_from_FASTA
#'
#' The function generates an FM_index structure starting
#'     from an input FASTA file.
#' The FM Index contains the suffix array (SA),
#'     the Burrows-Wheeler Transform (BWT),
#'     the Occurrencies matrix (Occ) and the Count array (C).
#'
#' The FASTA file is parsed using a Biostrings function (readDNAStringSet)
#'     in order to ensure a robust and reliable parsing method.
#'
#' The FM index is compressed by downsampling the suffix array,
#'     taking one every 32 suffixes.
#'
#' @param input path to the FASTA file
#' @param output path to the directory in which to save the data structures
#' @param save default = TRUE, can be modified
#'     to avoid saving data in the output directory
#' @param compress default = FALSE, can be modified
#'     to compress the FM index, by downsampling the suffix array
#' @return List-like object of type FM_index
#' @importFrom zoo na.locf
#' @importFrom Biostrings readDNAStringSet
#' @importFrom utils write.table
#' @examples
#' # Example creation of an FM index from a FASTA file
#' input_file <- system.file("extdata", "prova.txt", package = "FMIndexR")
#' output_path <- system.file("output", package = "FMIndexR")
#'
#' FM_index <- FM_index_from_FASTA(input_file,
#'                                 output_path,
#'                                 save = FALSE,
#'                                 compress = FALSE)
#' FM_index
#' @export
FM_index_from_FASTA <- function(input, output, save = TRUE, compress = FALSE){
    if (!file.exists(input) || (!dir.exists(output))) {
        stop("File or output directory NOT found!")
        } else {
        options(scipen = 999)
        fasta_data <- Biostrings::readDNAStringSet(input)
        if (length(fasta_data) > 1) {
            stop("multiFASTA files are not accepted!")
            } else {
                fasta_sequence <- as.character(fasta_data)
                fasta_sequence <- toupper(fasta_sequence)
                fasta_sequence_header <- names(fasta_data)
                fasta_sequence_header_path <- paste(output,"name.txt",
                                                    sep = "",collapse = "")
                message("FASTA sequence correctly read!")

                suffix_array <- .SuffixArray(fasta_sequence)
                BWT <- .BWTransform(suffix_array)
                if (compress) {
                    suffix_array <- suffix_array[seq(1,nrow(suffix_array),32),]}

                suffix_array_path <- paste(output,"suffix_array.txt",
                                            sep = "",collapse = "")
                BWT_path <- paste(output,"BWT.txt",
                                    sep = "",collapse = "")

                message("Suffix Array created!")
                message("BWT created!")

                occ_matrix <- .OccMatrix(BWT)
                occ_matrix_path <- paste(output,"occ_matrix.txt",
                                            sep = "",collapse = "")
                message("Occ matrix created!")

                c_array <- .CountArray(BWT)
                c_array_path <- paste(output,"c_array.txt",
                                        sep = "",collapse = "")
                message("C Array created!")

                if (save) {
                    .save(fasta_sequence_header,fasta_sequence_header_path,
                            suffix_array,suffix_array_path,BWT,BWT_path,
                            occ_matrix,occ_matrix_path,c_array,c_array_path)
                    message("Files correctly created!")}
                FM_index <- .FMIndex(fasta_sequence_header,suffix_array,
                                        BWT,occ_matrix,c_array)
                options(scipen = 0)
                return(FM_index)}}}
