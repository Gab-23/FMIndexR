#' FM_index_from_FASTA
#'
#' Generates FM_index structure starting from an input FASTA file
#'
#' FM index is composed of 4 different elements:
#'
#' - Suffix Array --> stored as a Dataframe, contains all suffixes and their indexes of the original string
#' - Burrows-Wheeler Transform --> a transformed string
#' - Occurrency matrix --> Dataframe storing the cumulative count table of characters in the BWT
#' - Count array --> Dataframe containing, for each character, the count of characters that are lexicographically smaller in the sorted BWT
#'
#' The FASTA file is parsed using a Biostrings function (readDNAStringSet) in order to ensure a robust and reliable parsing method
#' The FM Index returned by the function contains also the header of the FASTA file in order to keep informations related to the original sequence
#'
#' @param input path to the FASTA file
#' @param output path to the folder in which to save the data structures
#' @param save default = TRUE, can be modified to avoid saving data in the output path as separate files
#' @return Object of type FM_index (list-like object)
#' @importFrom zoo na.locf
#' @importFrom Biostrings readDNAStringSet
#' @importFrom utils write.table
#' @export
FM_index_from_FASTA <- function(input, output, save = TRUE){

  SuffixArray <- function(input_string) {

    if (typeof(input_string) != 'character') {
      stop('ERROR! Input sequence MUST be a string')

    } else if (nchar(input_string) == 0) {
      stop('ERROR! Empty sequence detected!')

    } else {
      special_char <- '$'
      complete_string <- paste(input_string,special_char, sep = '', collapse = '')
      complete_string_length <- nchar(complete_string)

      idx_array <- 0:(complete_string_length-1)

      suffix_array <- sapply(X = idx_array, function(i){
        actual_idx <- i
        used_idx <- i + 1
        suffix <- substr(complete_string, used_idx,complete_string_length)})

      suffix_df <- data.frame(idx = idx_array,
                              suffix = suffix_array)

      sorted_suffix_df <- suffix_df[order(suffix_array),]
      rownames(sorted_suffix_df) <- 0:(nrow(sorted_suffix_df)-1)
      class(sorted_suffix_df) <- c('SuffixArray','data.frame')
      return(sorted_suffix_df)
    }}
  BWTransform <- function(suffix_array) {

    logical <- any(class(suffix_array) != c('SuffixArray','data.frame'))

    if (logical) {
      stop('ERROR! Input MUST be a data.frame of class SuffixArray, use the dedicated function!')
    } else {
      original_sequence <- suffix_array[suffix_array$idx == 0,2]
      actual_idx <- suffix_array$idx
      used_idx <- ifelse(actual_idx == 0, nchar(original_sequence),actual_idx)

      BWT <- sapply(X = used_idx, function(i){
        idx <- i
        substr(original_sequence,idx,idx)})

      BWT <- paste(BWT, sep = '', collapse = '')
      class(BWT) <- c('BWT','character')
      return(BWT)
    }}
  OccMatrix <- function(bwt) {

    logical <- any(class(bwt) != c('BWT','character'))

    if (logical) {
      stop('ERROR! Input MUST be a string of class BWT, use the dedicated function!')
    } else {
      length_bwt <- nchar(bwt)
      bwt_letters <- strsplit(bwt, split = '')[[1]]
      unique_bwt_letters <- unique(bwt_letters)
      sorted_unique_bwt_letters <- sort(unique_bwt_letters)

      Occ_df <- data.frame(matrix(data = NA, nrow = length_bwt, ncol = length(sorted_unique_bwt_letters)))
      colnames(Occ_df) <- sorted_unique_bwt_letters

      for (char in colnames(Occ_df)) {
        idx_vector <- which(bwt_letters == char)
        Occ_df[idx_vector,char] <- 1:length(idx_vector)

        count_vec <- Occ_df[,char]
        count_vec <- zoo::na.locf(count_vec, na.rm = FALSE)
        count_vec[is.na(count_vec)] <- 0

        Occ_df[,char] <- count_vec

      }

      rownames(Occ_df) <- 0:(nrow(Occ_df)-1)
      class(Occ_df) <- c('Occ','data.frame')
      return(Occ_df)
    }}
  CountArray <- function(bwt) {

    logical <- any(class(bwt) != c('BWT','character'))

    if (logical) {
      stop('ERROR! Input MUST be a string of class BWT, use the dedicated function!')
    } else {
      bwt_char <- sort(unique((strsplit(bwt, split = '')[[1]])))
      sorted_bwt <- sort(strsplit(bwt, split = '')[[1]])

      count_df <- data.frame(matrix(data = NA, nrow = 1, ncol = length(bwt_char)))
      colnames(count_df) <- bwt_char

      for (char in bwt_char) {
        df_idx <- which(colnames(count_df) == char)
        char_idx <- which(sorted_bwt == char)[1] - 1
        count_df[1,df_idx] <- char_idx
      }

      rownames(count_df) <- 0
      class(count_df) <- c('C','data.frame')
      return(count_df)
    }}

  fasta_data <- tryCatch({Biostrings::readDNAStringSet(input)},
                         warning = function(w)
                         {stop('Non standard symbols detected in FASTA file!')})

  if (length(fasta_data) > 1) {
    stop('ERROR! multiFASTA files are not accepted!')
  } else {

    fasta_sequence <- as.character(fasta_data)
    fasta_sequence <- toupper(fasta_sequence)

    fasta_sequence_header <- names(fasta_data)
    fasta_sequence_header_path <- paste(output, 'name.txt', sep = '', collapse = '')
    print('FASTA sequence correctly read!')

    suffix_array <- SuffixArray(fasta_sequence)
    suffix_array_path <- paste(output, 'suffix_array.txt', sep = '', collapse = '')
    print('Suffix Array created!')

    BWT <- BWTransform(suffix_array)
    BWT_path <- paste(output, 'BWT.txt', sep = '', collapse = '')
    print('BWT created!')

    occ_matrix <- OccMatrix(BWT)
    occ_matrix_path <- paste(output, 'occ_matrix.txt', sep = '', collapse = '')
    print('Occ matrix created!')

    c_array <- CountArray(BWT)
    c_array_path <- paste(output, 'c_array.txt', sep = '', collapse = '')
    print('C Array created!')

    if (save) {
      writeLines(fasta_sequence_header, con = fasta_sequence_header_path)
      utils::write.table(suffix_array, file = suffix_array_path, sep = '\t', row.names = FALSE)
      writeLines(BWT, con = BWT_path)
      utils::write.table(occ_matrix, file = occ_matrix_path, sep = '\t', row.names = FALSE)
      utils::write.table(c_array, file = c_array_path, sep = '\t', row.names = FALSE)
    }

    FM_index <- list(fasta_sequence_header,suffix_array,BWT,occ_matrix,c_array)
    names(FM_index) <- c('SequenceName','SuffixArray', 'BWT', 'Occ', 'CountArray')
    class(FM_index) <- 'FM_index'

    return(FM_index)

  }}
