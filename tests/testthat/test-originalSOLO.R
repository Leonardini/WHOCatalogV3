test_that("originalSOLO classifies a variant as R when it is the sole variant in an R sample", {
  df  <- tibble(sample_id = 1L, variant = "V1", phenotype = "R", SOnly = FALSE)
  out <- originalSOLO(df, maxIter = 1L)
  expect_equal(out$class$mutation, "V1")
  expect_equal(out$class$class,    RCODE)
  expect_equal(out$class$Rcnt,     1L)
})

test_that("originalSOLO classifies a variant as S via SOnly flag at iteration 1", {
  df  <- tibble(sample_id = 1L, variant = "V1", phenotype = "S", SOnly = TRUE)
  out <- originalSOLO(df, maxIter = 1L)
  expect_equal(out$class$class, SCODE)
  expect_equal(out$class$iter,  1L)
})

test_that("originalSOLO leaves variants unclassified when no singletons exist within maxIter", {
  df <- tibble(
    sample_id = c(1L, 1L),
    variant   = c("V1", "V2"),
    phenotype = c("R", "R"),
    SOnly     = FALSE
  )
  out <- originalSOLO(df, maxIter = 1L, removeSOnly = FALSE)
  expect_true(all(out$class$class == UCODE))
})

test_that("originalSOLO with removeSOnly = FALSE skips the initial S-only classification", {
  df <- tibble(sample_id = 1L, variant = "V1", phenotype = "S", SOnly = TRUE)
  out_with    <- originalSOLO(df, maxIter = 1L, removeSOnly = TRUE)
  out_without <- originalSOLO(df, maxIter = 1L, removeSOnly = FALSE)
  expect_equal(out_with$class$class,    SCODE)
  expect_equal(out_without$class$class, SCODE)
})

test_that("originalSOLO with listIsolates = TRUE returns an isolate list", {
  df  <- tibble(sample_id = 1L, variant = "V1", phenotype = "R", SOnly = FALSE)
  out <- originalSOLO(df, maxIter = 1L, listIsolates = TRUE)
  expect_true("isolate" %in% names(out))
  expect_equal(out$isolate$variant,   "V1")
  expect_equal(out$isolate$sample_id, "1")
})

test_that("originalSOLO records Tcnt and Rcnt correctly for a solo R variant", {
  df <- tibble(
    sample_id = c(1L, 2L),
    variant   = c("V1", "V1"),
    phenotype = c("R", "S"),
    SOnly     = FALSE
  )
  out <- originalSOLO(df, maxIter = 1L)
  expect_equal(out$class$Tcnt, 2L)
  expect_equal(out$class$Rcnt, 1L)
})

test_that("originalSOLO uses the het column when present instead of SOnly", {
  df <- tibble(
    sample_id = 1L,
    variant   = "V1",
    phenotype = "S",
    het       = FALSE
  )
  out <- originalSOLO(df, maxIter = 1L)
  expect_equal(out$class$class, SCODE)
})
