library(testthat)
library(FMIndexR)

input_file <- system.file("extdata", "prova.txt", package = "FMIndexR")
output_path <- system.file("output", package = "FMIndexR")

FM_index <- suppressMessages(FM_index_from_FASTA(input_file, output_path, save = FALSE))
FM_index_comp <- suppressMessages(FM_index_from_FASTA(input_file, output_path, save = FALSE, compress = TRUE))

not_an_FM_index <- 2
not_a_pattern <- 2

example_search <- suppressMessages(BackwardSearch(FM_index, "CC", TRUE))
example_search_comp <- suppressMessages(BackwardSearch(FM_index_comp, "CC", TRUE))

example_search_2 <- suppressMessages(BackwardSearch(FM_index, "CCCCCCCCC", TRUE))

test_that("Error with no FM_index", {
  expect_error(BackwardSearch(not_an_FM_index, "AGA"))
})

test_that("Error with empty pattern", {
  expect_error(BackwardSearch(FM_index, ""))
})

test_that("Error with no pattern", {
  expect_error(BackwardSearch(FM_index, not_a_pattern))
})

test_that("Found pattern", {
  expect_equal(example_search[[2]], c(2,3,6))
})

test_that("Pattern not found", {
  expect_null(example_search_2)
})

test_that("Match compression", {
  expect_true(all(example_search$indexes==example_search_comp$indexes))
})

