make_sensspec_input <- function(n_samples = 5L, phenotype = "R") {
  tibble(
    sample_id      = seq_len(n_samples),
    drug           = "Isoniazid",
    variant        = "katG_p.Ser315Thr",
    gene           = "katG",
    mutation       = "p.Ser315Thr",
    phenotype      = phenotype,
    effect         = "missense_variant",
    tier           = 1L,
    `max(af)`      = 0.9,
    `max(quality)` = 2000L,
    pos1           = 315L,
    pos2           = NA_real_,
    position       = NA_character_,
    lineage        = "1"
  )
}

make_sensspec_catalogue <- function(td, grade = GRADES[1L]) {
  path <- file.path(td, "catalogue.csv")
  write_csv(
    tibble(drug = "Isoniazid", variant = "katG_p.Ser315Thr",
           Final_Confidence_Grading = grade),
    path
  )
  path
}

test_that("computeSensSpec returns a list with a fullTab element", {
  td <- withr::local_tempdir()
  result <- computeSensSpec(
    make_sensspec_input(),
    catalogueFile = make_sensspec_catalogue(td),
    safe = FALSE
  )
  expect_true("fullTab" %in% names(result))
})

test_that("computeSensSpec fullTab contains expected count columns", {
  td <- withr::local_tempdir()
  result <- computeSensSpec(
    make_sensspec_input(),
    catalogueFile = make_sensspec_catalogue(td),
    safe = FALSE
  )
  expect_true(all(c("drug", "TP", "TN", "FP", "FN") %in% colnames(result$fullTab)))
})

test_that("computeSensSpec places grade-1 R samples in Group_Regular 1", {
  td <- withr::local_tempdir()
  result <- computeSensSpec(
    make_sensspec_input(n_samples = 5L, phenotype = "R"),
    catalogueFile = make_sensspec_catalogue(td, grade = GRADES[1L]),
    safe = FALSE
  )
  inh <- result$fullTab[result$fullTab$drug == "Isoniazid" & grepl("Regular_1$", result$fullTab$group), ]
  expect_true(nrow(inh) > 0L)
  expect_true(all(inh$TP > 0L))
})
