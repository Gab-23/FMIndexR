library(testthat)
library(FMIndexR)

input_file_ok <- system.file("extdata", "prova.txt", package = "FMIndexR")
input_file_no_header <- system.file("extdata", "prova_no_header.txt", package = "FMIndexR")
input_file_warning <- system.file("extdata", "prova_warning.txt", package = "FMIndexR")
input_file_seq_vuota <- system.file("extdata", "prova_seq_vuota.txt", package = "FMIndexR")
input_file_multifasta <- system.file("extdata", "prova_multifasta.txt", package = "FMIndexR")

output_path <- system.file("output", package = "FMIndexR")

FM_index_ok <- suppressMessages(FM_index_from_FASTA(input_file_ok,
                                                    output_path,
                                                    save = FALSE))

test_that("Correct FM_index", {
  expect_s3_class(FM_index_ok, "FM_index")
})

test_that("Error with no header", {
  expect_error(FM_index_from_FASTA(input_file_no_header,
    output_path,
    save = FALSE
  ))
})

test_that("Error with seq vuota", {
  expect_error(suppressMessages(FM_index_from_FASTA(input_file_seq_vuota,
    output_path,
    save = FALSE
  )))
})

test_that("Warning with strange characters", {
  expect_warning(suppressMessages(FM_index_from_FASTA(input_file_warning,
    output_path,
    save = FALSE
  )))
})

test_that("Error with multifasta", {
  expect_error(FM_index_from_FASTA(input_file_multifasta,
                                   output_path,
                                   save = FALSE
  ))
})

test_that("File not found", {
  expect_error(FM_index_from_FASTA("not-a-file",
                                   output_path,
                                   save = FALSE
  ))
})

test_that("Dir not found", {
  expect_error(FM_index_from_FASTA(input_file_ok,
                                   "not-a-dir",
                                   save = FALSE
  ))
})
