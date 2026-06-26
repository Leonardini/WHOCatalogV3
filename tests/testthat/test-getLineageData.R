testthat::local_edition(3)

write_lineage_fixture <- function(tab) {
  dir <- withr::local_tempdir(.local_envir = parent.frame())
  countsDir <- file.path(dir, "extract", "lineages_counts")
  dir.create(countsDir, recursive = TRUE)
  write_csv(tab, file.path(countsDir, "lineages_counts.csv"))
  dir
}

make_lineage_row <- function(sample_id, lineage_numbering, final_af = 0.99) {
  tibble(sample_id = sample_id, position = sample_id, final_af = final_af,
         lineage_numbering = lineage_numbering)
}

test_that("getLineageData keeps lineages 1..LINEAGE_MAX and lumps the rest into Other", {
  fixture <- bind_rows(
    make_lineage_row(1, "lineage1.2.1"),
    make_lineage_row(2, "lineage4.1"),
    make_lineage_row(3, "lineage6"),
    make_lineage_row(4, "lineage7"),
    make_lineage_row(5, "lineage8.1"),
    make_lineage_row(6, "lineage9"),
    make_lineage_row(7, "lineageBOV"),
    make_lineage_row(8, "lineageLa1.2")
  )
  dir <- write_lineage_fixture(fixture)
  out <- getLineageData(dir, "extract", useSublineageData = FALSE)

  lineageOf <- function(id) out$lineage[out$sample_id == id]
  expect_equal(lineageOf(1), "1")
  expect_equal(lineageOf(2), "4")
  expect_equal(lineageOf(3), "6")
  expect_equal(lineageOf(4), LINEAGE_OTHER)
  expect_equal(lineageOf(5), LINEAGE_OTHER)
  expect_equal(lineageOf(6), LINEAGE_OTHER)
  expect_equal(lineageOf(7), LINEAGE_OTHER)
  expect_equal(lineageOf(8), LINEAGE_OTHER)
})

test_that("getLineageData preserves sub-lineage labels for in-range lineages", {
  fixture <- bind_rows(
    make_lineage_row(1, "lineage4.1.2"),
    make_lineage_row(2, "lineage2.2"),
    make_lineage_row(3, "lineage7.1")
  )
  dir <- write_lineage_fixture(fixture)
  out <- getLineageData(dir, "extract", useSublineageData = TRUE)

  lineageOf <- function(id) out$lineage[out$sample_id == id]
  expect_equal(lineageOf(1), "4.1")
  expect_equal(lineageOf(2), "2.2")
  expect_equal(lineageOf(3), LINEAGE_OTHER)
})
