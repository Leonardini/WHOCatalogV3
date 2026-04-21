test_that("Extracts positions correctly with valid input", {
  df <- tibble(mutation = c("A123T", "del10_ins-20", "ins5_7"))
  out <- extractPositions(df, "mutation", 2)
  
  expect_equal(out$pos1, c(123, 10, 5))
  expect_equal(out$pos2, c(NA_integer_, -20, 7))
})

test_that("Adds NA when fewer integers than K", {
  df <- tibble(mutation = c("A123T", "ins5"))
  out <- extractPositions(df, "mutation", 3)
  
  expect_equal(out$pos1, c(123, 5))
  expect_equal(out$pos2, c(NA_integer_, NA_integer_))
  expect_equal(out$pos3, c(NA_integer_, NA_integer_))
})

test_that("Throws error when a row has more than K integers", {
  df <- tibble(mutation = c("del5-10", "A1B2C3"))
  expect_error(
    extractPositions(df, "mutation", 2),
    "Rows 2 contain more than 2 integers"
  )
})

test_that("Handles NA and malformed strings", {
  df <- tibble(mutation = c(NA_integer_, "invalid", "A1"))
  out <- extractPositions(df, "mutation", 2)
  
  expect_equal(out$pos1, c(NA_integer_, NA_integer_, 1))
  expect_equal(out$pos2, c(NA_integer_, NA_integer_, NA_integer_))
})

test_that("Throws error for invalid inputs", {
  df <- tibble(mutation = "A123T")
  expect_error(extractPositions(df, "nonexistent", 2), "colName %in% names")
  expect_error(extractPositions(df, "mutation", -1), "maxNumber > 0")
})
