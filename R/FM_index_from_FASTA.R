#' FM_index_from_FASTA
#'
#' Generates FM_index structure starting from an input FASTA file
#'
#' The FASTA file is parsed using a Biostrings function (readDNAStringSet)
#' in order to ensure a robust and reliable parsing method
#'
#' @param input path to the FASTA file
#' @param output path to the folder in which to save the data structures
#' @param save default = TRUE, can be modified
#'     to avoid saving data in the output path as separate files
#' @return Object of type FM_index (list-like object)
#' @importFrom zoo na.locf
#' @importFrom Biostrings readDNAStringSet
#' @importFrom utils write.table
#' @examples
#' # example creation of an FM index from a FASTA file
#' input_file <- system.file("extdata", "prova.txt", package = "FMIndexR")
#' output_path <- system.file("output", package = "FMIndexR")
#' FM_index <- FM_index_from_FASTA(input_file, output_path, save = FALSE)
#' FM_index
#' @export
FM_index_from_FASTA <- function(input, output, save = TRUE) {
  SuffixArray <- function(input_string) {
    if (nchar(input_string) == 0) {
      stop("ERROR! Empty sequence detected!")
    } else {
      special_char <- "$"
      complete_string <- paste(input_string,
                               special_char,
                               sep = "", collapse = "")
      complete_string_length <- nchar(complete_string)

      idx_array <- 0:(complete_string_length - 1)

      suffix_array <- vapply(X = idx_array, function(i) {
        actual_idx <- i
        used_idx <- i + 1
        suffix <- substr(complete_string, used_idx, complete_string_length)
      }, character(1))

      suffix_df <- data.frame(
        idx = idx_array,
        suffix = suffix_array
      )

      sorted_suffix_df <- suffix_df[order(suffix_array), ]
      rownames(sorted_suffix_df) <- 0:(nrow(sorted_suffix_df) - 1)
      return(sorted_suffix_df)
    }
  }
  BWTransform <- function(suffix_array) {
    original_sequence <- suffix_array[suffix_array$idx == 0, 2]
    actual_idx <- suffix_array$idx
    used_idx <- ifelse(actual_idx == 0, nchar(original_sequence), actual_idx)
    BWT <- vapply(X = used_idx, function(i) {
      idx <- i
      substr(original_sequence, idx, idx)
    }, character(1))

    BWT <- paste(BWT, sep = "", collapse = "")
    return(BWT)
  }
  OccMatrix <- function(bwt) {
    length_bwt <- nchar(bwt)
    bwt_letters <- strsplit(bwt, split = "")[[1]]
    unique_bwt_letters <- unique(bwt_letters)
    sorted_unique_bwt_letters <- sort(unique_bwt_letters)

    Occ_df <- data.frame(matrix(data = NA,
                                nrow = length_bwt,
                                ncol = length(sorted_unique_bwt_letters)))
    colnames(Occ_df) <- sorted_unique_bwt_letters

    for (char in colnames(Occ_df)) {
      idx_vector <- which(bwt_letters == char)
      Occ_df[idx_vector, char] <- seq_along(idx_vector)

      count_vec <- Occ_df[, char]
      count_vec <- zoo::na.locf(count_vec, na.rm = FALSE)
      count_vec[is.na(count_vec)] <- 0

      Occ_df[, char] <- count_vec
    }

    rownames(Occ_df) <- 0:(nrow(Occ_df) - 1)
    return(Occ_df)
  }
  CountArray <- function(bwt) {
    bwt_char <- sort(unique((strsplit(bwt, split = "")[[1]])))
    sorted_bwt <- sort(strsplit(bwt, split = "")[[1]])

    count_df <- data.frame(matrix(data = NA, nrow = 1, ncol = length(bwt_char)))
    colnames(count_df) <- bwt_char

    for (char in bwt_char) {
      df_idx <- which(colnames(count_df) == char)
      char_idx <- which(sorted_bwt == char)[1] - 1
      count_df[1, df_idx] <- char_idx
    }

    rownames(count_df) <- 0
    return(count_df)
  }

  fasta_data <- Biostrings::readDNAStringSet(input)

  if (length(fasta_data) > 1) {
    stop("ERROR! multiFASTA files are not accepted!")
  } else {
    fasta_sequence <- as.character(fasta_data)
    fasta_sequence <- toupper(fasta_sequence)

    fasta_sequence_header <- names(fasta_data)
    fasta_sequence_header_path <- paste(output,
                                        "name.txt",
                                        sep = "", collapse = "")
    message("FASTA sequence correctly read!")

    suffix_array <- SuffixArray(fasta_sequence)
    suffix_array_path <- paste(output,
                               "suffix_array.txt",
                               sep = "", collapse = "")
    message("Suffix Array created!")

    BWT <- BWTransform(suffix_array)
    BWT_path <- paste(output, "BWT.txt", sep = "", collapse = "")
    message("BWT created!")

    occ_matrix <- OccMatrix(BWT)
    occ_matrix_path <- paste(output, "occ_matrix.txt", sep = "", collapse = "")
    message("Occ matrix created!")

    c_array <- CountArray(BWT)
    c_array_path <- paste(output, "c_array.txt", sep = "", collapse = "")
    message("C Array created!")

    if (save) {
      writeLines(fasta_sequence_header, con = fasta_sequence_header_path)
      utils::write.table(suffix_array,
                         file = suffix_array_path,
                         sep = "\t",
                         row.names = FALSE)
      writeLines(BWT, con = BWT_path)
      utils::write.table(occ_matrix,
                         file = occ_matrix_path,
                         sep = "\t",
                         row.names = FALSE)
      utils::write.table(c_array,
                         file = c_array_path,
                         sep = "\t",
                         row.names = FALSE)
    }

    FM_index <- list(fasta_sequence_header,
                     suffix_array,
                     BWT,
                     occ_matrix,
                     c_array)
    names(FM_index) <- c("SequenceName",
                         "SuffixArray",
                         "BWT",
                         "Occ",
                         "CountArray")
    class(FM_index) <- "FM_index"

    return(FM_index)
  }
}
