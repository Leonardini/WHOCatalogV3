testthat::local_edition(3)

test_that("recodes matched rows and leaves unmatched rows unchanged", {
  initTab  <- tibble(key = c("a", "b", "c"), val = c(1L, 2L, 3L))
  manualTab <- tibble(key = c("a", "c"), val = c(10L, 30L))
  out <- recodeValues(initTab, manualTab)
  expect_equal(nrow(out), 3L)
  expect_equal(out %>% arrange(key) %>% pull(val), c(10L, 2L, 30L))
})

test_that("returns full initTab when no keys match", {
  initTab  <- tibble(key = c("a", "b"), val = c(1L, 2L))
  manualTab <- tibble(key = "z", val = 99L)
  out <- recodeValues(initTab, manualTab)
  expect_equal(out, initTab)
})

test_that("errors when manualTab has fewer than 2 columns", {
  expect_snapshot(
    recodeValues(tibble(key = "a", val = 1), tibble(key = "a")),
    error = TRUE
  )
})

test_that("errors when manualTab has more than 2 columns", {
  expect_snapshot(
    recodeValues(tibble(key = "a", val = 1), tibble(a = 1, b = 2, c = 3)),
    error = TRUE
  )
})

test_that("errors when a manualTab column name ends with _first or _second", {
  expect_snapshot(
    recodeValues(tibble(key = "a", val = 1), tibble(key = "a", val_first = 10)),
    error = TRUE
  )
})
