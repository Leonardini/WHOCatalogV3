test_that("computePPVs computes PPV for each variant-drug pair", {
  df <- tibble(
    drug      = rep("DrugA", 3),
    variant   = c("V1", "V1", "V2"),
    het       = FALSE,
    phenotype = c("R", "S", "S"),
    sample_id = 1:3
  )
  out <- computePPVs(df, removeRM = FALSE, restrict = TRUE)
  expect_equal(nrow(out), 2L)
  expect_true("PPV" %in% colnames(out))
  expect_true("PPV_ub" %in% colnames(out))
  v1_ppv <- out$PPV[out$variant == "V1"]
  expect_equal(unname(v1_ppv), 0.5)
})

test_that("computePPVs with restrict = FALSE adds PPV columns to original table", {
  df <- tibble(
    drug      = "DrugA",
    variant   = "V1",
    het       = FALSE,
    phenotype = "R",
    sample_id = 1L
  )
  out <- computePPVs(df, removeRM = FALSE, restrict = FALSE)
  expect_equal(nrow(out), 1L)
  expect_true("PPV" %in% colnames(out))
})

test_that("computePPVs with removeRM = TRUE excludes samples marked as rm_in_sample", {
  df <- tibble(
    drug         = rep("DrugA", 2),
    variant      = "V1",
    het          = FALSE,
    phenotype    = c("R", "S"),
    sample_id    = 1:2,
    rm_in_sample = c(TRUE, FALSE)
  )
  out <- computePPVs(df, removeRM = TRUE, restrict = TRUE)
  expect_equal(unname(out$PPV[out$variant == "V1"]), 0.0)
})

test_that("computePPVs with solo = TRUE only counts single-variant samples", {
  df <- tibble(
    drug      = rep("DrugA", 4),
    variant   = c("V1", "V2", "V1", "V1"),
    het       = FALSE,
    phenotype = c("R", "R", "R", "S"),
    sample_id = c(1L, 1L, 2L, 3L)
  )
  out <- computePPVs(df, removeRM = FALSE, solo = TRUE, restrict = TRUE)
  v1_row <- out[out$variant == "V1", ]
  expect_equal(unname(v1_row$PPV), 0.5)
  expect_equal(nrow(out[out$variant == "V2", ]), 0L)
})

test_that("computePPVs treats het = TRUE rows as absent when computing totals", {
  df <- tibble(
    drug      = rep("DrugA", 2),
    variant   = "V1",
    het       = c(TRUE, FALSE),
    phenotype = c("R", "R"),
    sample_id = 1:2
  )
  out <- computePPVs(df, removeRM = FALSE, restrict = TRUE)
  expect_equal(unname(out$PPV), 1.0)
})
