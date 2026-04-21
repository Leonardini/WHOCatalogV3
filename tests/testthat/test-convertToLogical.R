test_that("converts NA to FALSE", {
  expect_false(convertToLogical(NA))
})

test_that("preserves TRUE", {
  expect_true(convertToLogical(TRUE))
})

test_that("preserves FALSE", {
  expect_false(convertToLogical(FALSE))
})

test_that("converts vector element-wise, replacing NA with FALSE", {
  expect_equal(convertToLogical(c(TRUE, NA, FALSE)), c(TRUE, FALSE, FALSE))
})

test_that("converts integer 1/0 to TRUE/FALSE", {
  expect_equal(convertToLogical(c(1L, 0L)), c(TRUE, FALSE))
})
