testthat::local_edition(3)

test_that("marks variants with MAF below minMAF as het when lowMAFHet = TRUE", {
  df <- tibble(`max(af)` = c(0.9, 0.5), variant = c("A", "B"))
  out <- applyThresholds(df, minMAF = 0.75, lowMAFHet = TRUE, minQ = NA)
  expect_equal(out$het, c(FALSE, TRUE))
})

test_that("marks variant named 'missing' as het regardless of MAF", {
  df <- tibble(`max(af)` = c(0.9, 0.9), variant = c("A", "missing"))
  out <- applyThresholds(df, minMAF = 0.75, lowMAFHet = TRUE, minQ = NA)
  expect_equal(out$het, c(FALSE, TRUE))
})

test_that("removes variants with MAF below minMAF when lowMAFHet = FALSE", {
  df <- tibble(`max(af)` = c(0.9, 0.5), variant = c("A", "B"))
  out <- applyThresholds(df, minMAF = 0.75, lowMAFHet = FALSE, minQ = NA)
  expect_equal(nrow(out), 1L)
  expect_equal(out$variant, "A")
})

test_that("marks variants with quality below minQ as het when lowQHet = TRUE", {
  df <- tibble(`max(af)` = c(0.9, 0.9), `max(quality)` = c(1200, 500), variant = c("A", "B"))
  out <- applyThresholds(df, minMAF = 0.75, lowMAFHet = TRUE, minQ = 1000, lowQHet = TRUE)
  expect_equal(out$het, c(FALSE, TRUE))
})

test_that("removes variants with quality below minQ when lowQHet = FALSE", {
  df <- tibble(`max(af)` = c(0.9, 0.9), `max(quality)` = c(1200, 500), variant = c("A", "B"))
  out <- applyThresholds(df, minMAF = NA, minQ = 1000, lowQHet = FALSE)
  expect_equal(nrow(out), 1L)
  expect_equal(out$variant, "A")
})

test_that("skips MAF filtering and prints warning when max(af) column is absent", {
  df <- tibble(variant = c("A", "B"))
  expect_snapshot(applyThresholds(df, minMAF = 0.75, minQ = NA))
})

test_that("skips quality filtering and prints warning when max(quality) column is absent", {
  df <- tibble(`max(af)` = c(0.9, 0.9), variant = c("A", "B"))
  expect_snapshot(applyThresholds(df, minMAF = NA, minQ = 1000))
})
