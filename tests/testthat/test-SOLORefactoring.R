testthat::local_edition(3)

test_that("SOLO classifies variants correctly using SOnly markers", {
  inputTab <- tibble(
    drug      = rep("Dummy", 7),
    sample_id = rep(1:4, c(3, 1, 2, 1)),
    variant   = LETTERS[c(1, 2, 3, 1, 3, 4, 2)],
    phenotype = rep(c("R", "S"), c(4, 3)),
    SOnly     = c(rep(FALSE, 5), TRUE, FALSE)
  )
  result <- runSOLOPipeline(inputTab, 2L, listIsolates = TRUE)
  expect_equal(result$classes$class, c("R", "S", "S", "S"))
  expect_equal(result$classes$iter,  c(1L, 1L, 2L, 1L))
  expect_equal(result$classes$Scnt,  c(0L, 1L, 1L, 0L))
  expect_equal(result$classes$Rcnt,  c(1L, 0L, 0L, 0L))
  expect_equal(result$isolates$variant,   LETTERS[1:3])
  expect_equal(result$isolates$sample_id, c(2L, 4L, 3L))
})

test_that("SOLO classifies variants correctly using the het column instead of SOnly", {
  inputTab <- tibble(
    drug      = rep("Dummy", 7),
    sample_id = rep(1:4, c(3, 1, 2, 1)),
    variant   = LETTERS[c(1, 2, 3, 1, 3, 4, 2)],
    phenotype = rep(c("R", "S"), c(4, 3)),
    SOnly     = FALSE
  )
  result <- runSOLOPipeline(inputTab, 2L, listIsolates = TRUE)
  expect_equal(result$classes$class, c("R", "S", "U", "U"))
  expect_equal(result$classes$iter,  c(1L, 1L, 0L, 0L))
  expect_equal(result$classes$Scnt,  c(0L, 1L, 0L, 0L))
  expect_equal(result$classes$Rcnt,  c(1L, 0L, 0L, 0L))
  expect_equal(result$isolates$variant,   LETTERS[1:2])
  expect_equal(result$isolates$sample_id, c(2L, 4L))
})

test_that("SOLO handles a tricky case with mixed phenotypes and a heterozygous call", {
  inputTab <- tibble(
    drug      = rep("Dummy", 6),
    sample_id = rep(1:3, each = 2),
    variant   = rep(c("x", "y"), 3),
    phenotype = c(rep("S", 4), rep("R", 2)),
    het       = c(rep(FALSE, 5), TRUE)
  )
  result <- runSOLOPipeline(inputTab, 2L, listIsolates = TRUE)
  expect_equal(result$classes$class, c("R", "S"))
  expect_equal(result$classes$iter,  c(2L, 1L))
  expect_equal(result$classes$Scnt,  c(2L, 0L))
  expect_equal(result$classes$Rcnt,  c(1L, 0L))
  expect_equal(result$isolates$variant,   rep("x", 3))
  expect_equal(result$isolates$sample_id, 1:3)
})

test_that("runSOLOPipeline errors when required column names are wrong", {
  expect_snapshot(
    runSOLOPipeline(tibble(a = 1, b = 2, c = 3, d = 4), maxIter = 1L),
    error = TRUE
  )
})