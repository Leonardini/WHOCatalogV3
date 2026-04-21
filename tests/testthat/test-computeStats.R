test_that("safeBinomTest returns list with NA fields when y = 0", {
  res <- safeBinomTest(0, 0)
  expect_true(is.na(res$estimate))
  expect_true(all(is.na(res$conf.int)))
})

test_that("safeBinomTest returns binom.test result when y > 0", {
  res <- safeBinomTest(3, 10)
  expect_false(is.na(res$estimate))
  expect_equal(unname(res$estimate), 0.3)
})

test_that("safeFisherTest returns list with NA p.value when first-row margin is zero", {
  expect_true(is.na(safeFisherTest(0, 0, 5, 5)$p.value))
})

test_that("safeFisherTest returns list with NA p.value when first-column margin is zero", {
  expect_true(is.na(safeFisherTest(0, 5, 0, 5)$p.value))
})

test_that("safeFisherTest returns list with NA p.value when second-row margin is zero", {
  expect_true(is.na(safeFisherTest(5, 5, 0, 0)$p.value))
})

test_that("safeFisherTest returns list with NA p.value when second-column margin is zero", {
  expect_true(is.na(safeFisherTest(5, 0, 5, 0)$p.value))
})

test_that("safeFisherTest returns fisher.test result for a valid non-degenerate table", {
  res <- safeFisherTest(10, 2, 3, 15)
  expect_false(is.na(res$p.value))
  expect_lt(res$p.value, 0.05)
})

test_that("quietMax suppresses warnings and returns -Inf for all-NA input", {
  expect_equal(quietMax(c(NA_real_, NA_real_), na.rm = TRUE), -Inf)
})

test_that("quietMax returns the correct maximum for normal numeric input", {
  expect_equal(quietMax(c(1, 3, 2)), 3)
})
