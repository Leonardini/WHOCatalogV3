testthat::local_edition(3)

make_ac_row <- function(drug, gene, mutation, effect = "missense_variant") {
  tibble(
    drug     = drug,
    gene     = gene,
    variant  = paste0(gene, "_", mutation),
    mutation = mutation,
    effect   = effect
  )
}

test_that("applyCatalogue errors when required columns are missing", {
  tf <- withr::local_tempfile(fileext = ".csv")
  write_csv(tibble(drug = "x", variant = "y", Final_Confidence_Grading = GRADES[3L]), tf)
  expect_snapshot(
    applyCatalogue(tibble(drug = "Isoniazid"), tf, minMAF = NA, minQ = NA),
    error = TRUE
  )
})

test_that("applyCatalogue assigns the catalogue grade to a matching variant", {
  tf <- withr::local_tempfile(fileext = ".csv")
  write_csv(
    tibble(drug = "Isoniazid", variant = "katG_p.Ser315Thr",
           Final_Confidence_Grading = GRADES[1L]),
    tf
  )
  expect_snapshot(
    applyCatalogue(make_ac_row("Isoniazid", "katG", "p.Ser315Thr"),
                   tf, minMAF = NA, minQ = NA, LoF = FALSE)
  )
})

test_that("applyCatalogue leaves Final as NA for variants not in the catalogue", {
  tf <- withr::local_tempfile(fileext = ".csv")
  write_csv(
    tibble(drug = "Isoniazid", variant = "katG_p.Ser315Thr",
           Final_Confidence_Grading = GRADES[1L]),
    tf
  )
  expect_snapshot(
    applyCatalogue(make_ac_row("Isoniazid", "katG", "p.Arg463Leu"),
                   tf, minMAF = NA, minQ = NA, LoF = FALSE)
  )
})

test_that("applyCatalogue sets Final = 2 for unmatched RRDR non-silent variants", {
  tf <- withr::local_tempfile(fileext = ".csv")
  write_csv(
    tibble(drug = "Rifampicin", variant = "rpoB_p.Leu430Pro",
           Final_Confidence_Grading = GRADES[1L]),
    tf
  )
  expect_snapshot(
    applyCatalogue(make_ac_row("Rifampicin", "rpoB", "p.Ser450Leu"),
                   tf, minMAF = NA, minQ = NA, LoF = FALSE)
  )
})

test_that("applyCatalogue sets Final = 2 for unmatched LoF variants when pooled LoF is graded 1 or 2", {
  tf <- withr::local_tempfile(fileext = ".csv")
  write_csv(
    tibble(drug = "Isoniazid", variant = "katG_LoF",
           Final_Confidence_Grading = GRADES[1L]),
    tf
  )
  expect_snapshot(
    applyCatalogue(make_ac_row("Isoniazid", "katG", "p.Ala12Val", effect = "frameshift"),
                   tf, minMAF = NA, minQ = NA, LoF = TRUE)
  )
})

test_that("applyCatalogue does not apply the LoF rule when LoF = FALSE", {
  tf <- withr::local_tempfile(fileext = ".csv")
  write_csv(
    tibble(drug = "Isoniazid", variant = "katG_LoF",
           Final_Confidence_Grading = GRADES[1L]),
    tf
  )
  expect_snapshot(
    applyCatalogue(make_ac_row("Isoniazid", "katG", "p.Ala12Val", effect = "frameshift"),
                   tf, minMAF = NA, minQ = NA, LoF = FALSE)
  )
})
