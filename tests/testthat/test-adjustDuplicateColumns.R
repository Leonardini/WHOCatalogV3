testthat::local_edition(3)

test_that("takes .x value when .y is NA, and .y value when .x is NA", {
  df <- tibble(col.x = c(1, 2, NA), col.y = c(NA, NA, 3))
  out <- adjustDuplicateColumns(df, warn = FALSE)
  expect_equal(out$col, c(1, 2, 3))
  expect_false("col.x" %in% colnames(out))
  expect_false("col.y" %in% colnames(out))
})

test_that("handles multiple duplicate column pairs", {
  df <- tibble(a.x = c(1, NA), a.y = c(NA, 2), b.x = c("p", "q"), b.y = c(NA, NA))
  out <- adjustDuplicateColumns(df, warn = FALSE)
  expect_equal(out$a, c(1, 2))
  expect_equal(out$b, c("p", "q"))
})

test_that("sums numeric columns when add = TRUE", {
  df <- tibble(col.x = c(1, 2, NA), col.y = c(3, NA, 4))
  out <- adjustDuplicateColumns(df, warn = FALSE, add = TRUE)
  expect_equal(out$col, c(4, 2, 4))
})

test_that("concatenates character columns with '; ' when add = TRUE and both non-NA", {
  df <- tibble(col.x = c("a", NA, "c"), col.y = c("b", "d", NA))
  out <- adjustDuplicateColumns(df, warn = FALSE, add = TRUE)
  expect_equal(out$col, c("a; b", "d", "c"))
})

test_that("prints warning for conflicting non-NA entries when warn = TRUE", {
  df <- tibble(col.x = c(1, 2), col.y = c(3, 2))
  expect_snapshot(adjustDuplicateColumns(df, warn = TRUE))
})
