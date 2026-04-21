test_that("All values are consistent", {
  df <- data.frame(
    ID = c(1,1,2,2),
    group = c("A", "A", "B", "B"),
    val1 = c("x", "x", "y", "y"),
    val2 = c(10,10,20,20)
  )
  res <- testConsistent(df, groupingVars = "group", consistentVars = c("val1", "val2"))
  expect_true(all(res$checks))
  expect_true(all(sapply(res$problems, is.null)))
})

test_that("Detects inconsistencies correctly", {
  df <- data.frame(
    ID = c(1,1,2,2),
    group = c("A", "A", "B", "B"),
    val1 = c("x", "x", "y", "z"),
    val2 = c(10,11,20,20)
  )
  res <- testConsistent(df, groupingVars = "group", consistentVars = c("val1", "val2"))
  expect_false(all(res$checks))
  expect_false(is.null(res$problems$val1))
  expect_false(is.null(res$problems$val2))
})

test_that("Returns TRUE for empty consistentVars", {
  df <- data.frame(a = 1:3, b = letters[1:3])
  expect_true(testConsistent(df, groupingVars = "a", consistentVars = character(0)))
})

test_that("Stops on invalid column names", {
  df <- data.frame(a = 1:3, b = letters[1:3])
  expect_error(testConsistent(df, groupingVars = "z", consistentVars = "b"))
  expect_error(testConsistent(df, groupingVars = "a", consistentVars = "z"))
})
