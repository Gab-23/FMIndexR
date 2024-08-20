#' BackwardSearch
#'
#' The function looks for patterns inside the original sequence,
#'     using its FM_index and performing a
#'     backward search procedure.
#'
#' If the Suffix Array is compressed Last-to-First mapping
#'     will be applied to reconstruct missing values
#'
#' Querying a sequence using a compressed FM index may take
#'    longer than using the uncompressed one
#'
#' When store_elems = TRUE the original sequence has to be retrieved,
#'    if not stored in the suffix array it has to be rebuilt from the BWT
#'    this operation may require some time for longer sequences
#'
#' @param FM_index object of class FM_index obtained using FM_index_from_FASTA
#' @param pattern non-empty string containing the pattern to look for
#' @param store_elems default = FALSE, can be modified to save
#'     the original sequence,
#'     the pattern indexes and the pattern string
#'     in case a manual check wants to be made
#' @return Shows a report of the number of patterns found and their indexes
#'     in the original sequence. The indexes are also returned.
#'     If store_elems = TRUE a more detailed report for manual check
#'     can be returned.
#'     If no pattern is found, NULL is returned
#' @importFrom IRanges reverse
#' @importFrom methods is
#' @examples
#' # Example pattern search using BackwardSearch and an FM Index
#' input_file <- system.file("extdata", "prova.txt", package = "FMIndexR")
#' output_path <- system.file("output", package = "FMIndexR")
#'
#' FM_index <- FM_index_from_FASTA(input_file,
#'                                 output_path,
#'                                 save = FALSE,
#'                                 compress = FALSE)
#'
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

                BWT_array <- unique(strsplit(BWT,split = "")[[1]])
                BWT_array <- BWT_array[BWT_array != "$"]

                pattern <- toupper(pattern)
                reversed_pattern <- IRanges::reverse(pattern)
                reversed_pattern_array <- strsplit(reversed_pattern,
                                                    split = "")[[1]]

                logical_2 <- all(unique(reversed_pattern_array) %in% BWT_array)

                if (!logical_2) {
                    message("Pattern NOT found")
                    return(NULL)}

                range_and_count <- .getRangeandCount(reversed_pattern_array,
                                                        Occ,C,BWT)
                if (range_and_count$num_pattern == 0) {
                    message("Pattern NOT found")
                    return(NULL)
                    } else {
                        final_range <- range_and_count$range
                        num_pattern <- range_and_count$num_pattern}

                indexes <- .getIndexes(BWT,SA,C,Occ,final_range,num_pattern)
                indexes <- sort(unlist(indexes))
                indexes_str <- paste(indexes, collapse = ", ")

                message(sprintf('Pattern: %s',pattern))
                message(sprintf("%d pattern(s) found",num_pattern))
                message(sprintf("Index(es): %s ",indexes_str))

                if (store_elems) {
                    original_sequence <- .get_original_sequence(SA,BWT,Occ,C)
                    return(list(sequence = original_sequence,
                                indexes = indexes,pattern = pattern))}
                    else {return(indexes)}
                }}
