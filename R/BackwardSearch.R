#' BackwardSearch
#'
#' Looks for patterns inside the original sequence, using its FM_index and performing a backward search procedure for finding patterns
#'
#' @param FM_index object of class FM_index obtained using FM_index_from_FASTA
#' @param pattern string containing the pattern to look for
#' @param store_elems default = FALSE, can be modified to save the original sequence, the pattern indexes and the pattern in case a manual check wants to be made
#' @return Prints a report of the number of patterns and their location, if store_elems = TRUE a report for manual check can be returned
#' @importFrom IRanges reverse
#' @importFrom methods is
#' @examples
#' # example creation of pattern search using BackwardSearch and an FM Index
#' input_file <- system.file("extdata", "prova.txt", package = "FMIndexR")
#' output_path <- system.file("output", package = "FMIndexR")
#' FM_index <- FM_index_from_FASTA(input_file, output_path, save = FALSE)
#' BackwardSearch(FM_index,"ATC")
#' @export
BackwardSearch <- function(FM_index, pattern, store_elems = FALSE) {

  if (!is(pattern,'character') || (nchar(pattern) == 0)) {
    stop('ERROR! Pattern MUST be a non-empty string')
  } else {
    if (!is(FM_index, "FM_index")) {
      stop('ERROR! FM index MUST be of class FM_index!')
    } else {

      SA <- FM_index$SuffixArray
      BWT <- FM_index$BWT
      Occ <- FM_index$Occ
      C <- FM_index$CountArray

      SA$no_money <-substr(SA$suffix, 1, nchar(SA$suffix) - 1)

      original_sequence <- SA[SA$idx == 0,]$suffix
      original_sequence <- substr(original_sequence, 1, nchar(original_sequence) - 1)
      original_sequence_array <- strsplit(original_sequence, split = '')[[1]]

      pattern <- toupper(pattern)
      reversed_pattern <- IRanges::reverse(pattern)
      reversed_pattern_array <- strsplit(reversed_pattern, split = '')[[1]]

      logical_2 <- length(setdiff(reversed_pattern_array,original_sequence_array)) > 0

      if (logical_2) {
        return('Pattern NOT found')
      }

      if (any(SA$no_money == pattern)) {
        match <- which(SA$no_money == pattern)
        match <- rownames(SA)[match]
      } else {match <- ''}

      start <- 1
      end <- (nchar(BWT)-1)

      for (char in reversed_pattern_array) {

        start <- C[as.character(0),char] + Occ[as.character(start-1),char]
        end <- C[as.character(0),char] + Occ[as.character(end),char] -1

        if ((start > end) & (nchar(match) == 0)) {
          return('Pattern NOT found')
          break
        }
      }

      if (nchar(match) != 0) {
        start <- as.integer(match)}

      final_range <- as.character((start):(end))

      num_pattern <- length(final_range)
      indexes <- sort(as.vector(SA[final_range,]$idx))
      indexes_str <- paste(indexes, collapse = ", ")

      if (nchar(original_sequence) <= 3000) {
        print(paste('Original Sequence: ', original_sequence, sep = ''))
      }

      print(paste(num_pattern, ' pattern(s) found', sep = ''))
      print(paste('Index(es): ', indexes_str, sep = ''))

      if (store_elems) {
        return(list(original_sequence,indexes,pattern))
      }

    }}}
