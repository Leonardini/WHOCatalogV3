test_that("includes sample-drug pairs where all phenotypes are S", {
  df <- tibble(
    sample_id = c(1L, 1L),
    drug      = c("A", "A"),
    phenotype = c("S", "S"),
    class     = c(NA_character_, NA_character_)
  )
  out <- findUnexplained(df)
  expect_equal(nrow(out), 2L)
})

test_that("includes R-phenotype sample-drug pairs where all classes are NA or S", {
  df <- tibble(
    sample_id = c(1L, 1L),
    drug      = c("A", "A"),
    phenotype = c("R", "R"),
    class     = c("S", NA_character_)
  )
  out <- findUnexplained(df)
  expect_equal(nrow(out), 2L)
})

test_that("excludes sample-drug pairs where any class is R", {
  df <- tibble(
    sample_id = c(1L, 1L),
    drug      = c("A", "A"),
    phenotype = c("R", "R"),
    class     = c("R", NA_character_)
  )
  out <- findUnexplained(df)
  expect_equal(nrow(out), 0L)
})

test_that("handles mixed sample-drug pairs correctly, keeping only unexplained ones", {
  df <- tibble(
    sample_id = c(1L, 2L),
    drug      = c("A", "A"),
    phenotype = c("S", "R"),
    class     = c(NA_character_, "R")
  )
  out <- findUnexplained(df)
  expect_equal(nrow(out), 1L)
  expect_equal(out$sample_id, 1L)
})

test_that("does not retain the satisfying column in the output", {
  df <- tibble(sample_id = 1L, drug = "A", phenotype = "S", class = NA_character_)
  out <- findUnexplained(df)
  expect_false("satisfying" %in% colnames(out))
})
