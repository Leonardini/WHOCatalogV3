test_that("returns dominant broad lineage per sample", {
  df <- tibble(
    sample_id        = c(1, 1, 1, 2, 2),
    position         = 1:5,
    lineage_numbering = c("lineage4.1.2", "lineage4.2", "lineage2.2", "lineage4.1", "lineage2.1"),
    final_af         = c(0.9, 0.8, 0.8, 0.9, 0.9)
  )
  out <- computeDominantLineage(df, mafThresholds = 0.75, subLineage = FALSE)[["maf0.75"]]

  expect_equal(nrow(out), 2)
  s1 <- out[out$sample_id == 1, ]
  expect_equal(s1$dominant_lineage, "4")
  expect_equal(s1$ties, "")
  s2 <- out[out$sample_id == 2, ]
  expect_equal(s2$dominant_lineage, "2")
  expect_equal(s2$ties, "4")
})

test_that("returns dominant sub-lineage per sample", {
  df <- tibble(
    sample_id        = c(1, 1, 1, 2),
    position         = 1:4,
    lineage_numbering = c("lineage4.1.2", "lineage4.1.3", "lineage2.2", "lineage4.1"),
    final_af         = c(0.9, 0.9, 0.9, 0.9)
  )
  out <- computeDominantLineage(df, mafThresholds = 0.75, subLineage = TRUE)[["maf0.75"]]

  s1 <- out[out$sample_id == 1, ]
  expect_equal(s1$dominant_lineage, "4.1")
  expect_equal(s1$ties, "")
  s2 <- out[out$sample_id == 2, ]
  expect_equal(s2$dominant_lineage, "4.1")
  expect_equal(s2$ties, "")
})

test_that("samples below MAF threshold appear with empty lineage", {
  df <- tibble(
    sample_id        = c(1, 2),
    position         = 1:2,
    lineage_numbering = c("lineage4.1", "lineage2.1"),
    final_af         = c(0.5, 0.9)
  )
  out <- computeDominantLineage(df, mafThresholds = 0.75, subLineage = FALSE)[["maf0.75"]]

  s1 <- out[out$sample_id == 1, ]
  expect_equal(s1$dominant_lineage, "")
  expect_equal(s1$ties, "")
  s2 <- out[out$sample_id == 2, ]
  expect_equal(s2$dominant_lineage, "2")
})

test_that("errors on missing required columns", {
  df <- tibble(sample_id = 1, lineage_numbering = "lineage4")
  expect_error(computeDominantLineage(df), "Missing required columns: final_af")
})
