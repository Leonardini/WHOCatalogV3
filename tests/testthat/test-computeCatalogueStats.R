testthat::local_edition(3)

make_stats_df <- function() {
  tibble(
    drug         = "DrugA",
    variant      = "V1",
    present_R    = 10L,
    present_S    = 2L,
    absent_R     = 3L,
    absent_S     = 15L,
    SOLO_R       = 8L,
    SOLO_S       = 1L,
    correctAll   = TRUE,
    correctSOLO  = TRUE
  )
}

test_that("computeCatalogueStats prints the number of rows being processed", {
  expect_snapshot(computeCatalogueStats(make_stats_df()))
})

test_that("computeCatalogueStats returns one row with key stat columns", {
  out <- capture.output(result <- computeCatalogueStats(make_stats_df()))
  expect_equal(nrow(result), 1L)
  expect_true(all(c("PPV", "PPV_SOLO", "PPVc_SOLO", "Sens", "Spec", "OR", "OR_SOLO") %in% colnames(result)))
  expect_true(all(c("OR_pval_FDR_sig", "OR_SOLO_pval_FDR_sig") %in% colnames(result)))
})

test_that("computeCatalogueStats removes rows where SOLO and absent counts are all zero", {
  df <- tibble(
    drug        = "DrugA",
    variant     = c("V1", "V2"),
    present_R   = c(10L, 5L),
    present_S   = c(2L,  3L),
    absent_R    = c(3L,  0L),
    absent_S    = c(15L, 0L),
    SOLO_R      = c(8L,  0L),
    SOLO_S      = c(1L,  0L),
    correctAll  = TRUE,
    correctSOLO = TRUE
  )
  out <- capture.output(result <- computeCatalogueStats(df))
  expect_equal(nrow(result), 1L)
  expect_equal(result$variant, "V1")
})

test_that("computeCatalogueStats marks significant OR_SOLO when evidence is strong", {
  out <- capture.output(result <- computeCatalogueStats(make_stats_df()))
  expect_true(result$OR_SOLO_pval_FDR_sig)
})
