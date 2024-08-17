#' BackwardSearch
#'
#' The function looks for patterns inside the original sequence,
#'     using its FM_index and performing a
#'     backward search procedure.
#'
#' If the Suffix Array is compressed Last-to-First mapping
#'     will be applied to reconstruct missing values
#'
#' Please note that querying a sequence using a compressed
#'     FM index will take longer than using the uncompressed one
#'
#' @param FM_index object of class FM_index obtained using FM_index_from_FASTA
#' @param pattern non-empty string containing the pattern to look for
#' @param store_elems default = FALSE, can be modified to save
#'     the original sequence,
#'     the pattern indexes and the pattern string
#'     in case a manual check wants to be made
#' @return Shows a report of the number of patterns found and their location.
#'     If store_elems = TRUE a report for manual check can be returned.
#'     If no pattern is found, NULL is returned
#' @importFrom IRanges reverse
#' @importFrom methods is
#' @examples
#' # Example pattern search using BackwardSearch and an FM Index
#' input_file <- system.file("extdata", "prova.txt", package = "FMIndexR")
#' output_path <- system.file("output", package = "FMIndexR")
#' FM_index <- FM_index_from_FASTA(input_file, output_path, save = FALSE)
#' BackwardSearch(FM_index, "CC")
#' @export
BackwardSearch <- function(FM_index, pattern, store_elems = FALSE) {
    if (!is(pattern, "character") || (nchar(pattern) == 0)) {
        stop("Pattern MUST be a non-empty string")
        } else if (!is(FM_index, "FM_index")) {
            stop("FM index MUST be of class FM_index!")
            } else {

                options(scipen = 999)
                SA <- FM_index$SuffixArray
                BWT <- FM_index$BWT
                Occ <- FM_index$Occ
                C <- FM_index$CountArray

                original_sequence <- .reverseBWT(BWT,C,Occ)
                original_sequence_array <- strsplit(original_sequence,
                                                    split = "")[[1]]
                pattern <- toupper(pattern)
                reversed_pattern <- IRanges::reverse(pattern)
                reversed_pattern_array <- strsplit(reversed_pattern,
                                                    split = "")[[1]]

                logical_2 <- length(setdiff(reversed_pattern_array,
                                            original_sequence_array)) > 0
                if (logical_2) {
                    message("Pattern NOT found")
                    return(NULL)}

                range <- .getRange(reversed_pattern_array,Occ,C,BWT)
                start <- range$start
                end <- range$end

                if (start == '/' && end == '/') {
                    message("Pattern NOT found")
                    return(NULL)}

                final_range <- as.character((start):(end))
                num_pattern <- length(final_range)
                indexes <- .getIndexes(BWT,SA,C,Occ,final_range,num_pattern)
                indexes <- sort(unlist(indexes))
                indexes_str <- paste(indexes, collapse = ", ")

                message(sprintf('Pattern: %s',pattern))
                message(sprintf("%d pattern(s) found",num_pattern))
                message(sprintf("Index(es): %s ",indexes_str))

                if (store_elems) {
                    return(list(sequence = original_sequence,
                                indexes = indexes,pattern = pattern))}}}
