test_that("removes all rows for a sample-drug pair that has at least one het", {
  df <- tibble(
    sample_id = c(1L, 1L, 2L, 2L),
    drug      = rep("A", 4),
    het       = c(FALSE, TRUE, FALSE, FALSE)
  )
  out <- removeSamplesWithHets(df)
  expect_equal(nrow(out), 2L)
  expect_equal(unique(out$sample_id), 2L)
})

test_that("keeps all rows when no sample-drug pair contains a het", {
  df <- tibble(
    sample_id = c(1L, 2L),
    drug      = rep("A", 2),
    het       = c(FALSE, FALSE)
  )
  out <- removeSamplesWithHets(df)
  expect_equal(nrow(out), 2L)
})

test_that("treats het status independently per sample-drug combination", {
  df <- tibble(
    sample_id = c(1L, 1L, 2L, 2L),
    drug      = c("A", "A", "A", "B"),
    het       = c(TRUE, FALSE, FALSE, FALSE)
  )
  out <- removeSamplesWithHets(df)
  expect_equal(nrow(out), 2L)
  expect_equal(unique(out$sample_id), 2L)
})

test_that("does not retain the satisfying column in the output", {
  df <- tibble(sample_id = 1L, drug = "A", het = FALSE)
  out <- removeSamplesWithHets(df)
  expect_false("satisfying" %in% colnames(out))
})
