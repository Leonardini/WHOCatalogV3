make_final_catalog <- function() {
  tibble(
    drug                              = c("Isoniazid", "Rifampicin"),
    variant                           = c("katG_p.Ser315Thr", "rpoB_p.Ser450Leu"),
    Initial_Confidence_Grading        = c(1L, 1L),
    Supplementary_Grading_Considerations = c("Rule: LoF", "Rule: LoF"),
    Final_Confidence_Grading          = c(1L, 1L)
  )
}

make_manual_checks <- function() {
  tibble(
    drug                              = "Isoniazid",
    variant                           = "katG_p.Ser315Thr",
    Supplementary_Grading_Considerations = "Evidence: manual",
    Final_Confidence_Grading          = 2L
  )
}

test_that("applyManualChecks overrides Final_Confidence_Grading for matched variants", {
  result <- applyManualChecks(make_final_catalog(), make_manual_checks())
  inh_row <- result[result$drug == "Isoniazid", ]
  expect_equal(inh_row$Final_Confidence_Grading, 2L)
})

test_that("applyManualChecks overrides Supplementary_Grading_Considerations for matched variants", {
  result <- applyManualChecks(make_final_catalog(), make_manual_checks())
  inh_row <- result[result$drug == "Isoniazid", ]
  expect_equal(inh_row$Supplementary_Grading_Considerations, "Evidence: manual")
})

test_that("applyManualChecks leaves unmatched variants unchanged", {
  result <- applyManualChecks(make_final_catalog(), make_manual_checks())
  rif_row <- result[result$drug == "Rifampicin", ]
  expect_equal(rif_row$Final_Confidence_Grading, 1L)
  expect_equal(rif_row$Supplementary_Grading_Considerations, "Rule: LoF")
})

test_that("applyManualChecks adds rows present only in manualCheckResults", {
  manual <- tibble(
    drug                              = "Pyrazinamide",
    variant                           = "pncA_p.Trp119Stop",
    Supplementary_Grading_Considerations = "Evidence: manual",
    Final_Confidence_Grading          = 1L
  )
  result <- applyManualChecks(make_final_catalog(), manual)
  expect_true("pncA_p.Trp119Stop" %in% result$variant)
})

test_that("applyManualChecks with empty manualCheckResults returns catalog unchanged", {
  empty <- tibble(drug = character(), variant = character(),
                  Supplementary_Grading_Considerations = character(),
                  Final_Confidence_Grading = integer())
  result <- applyManualChecks(make_final_catalog(), empty)
  expect_equal(nrow(result), 2L)
  expect_equal(result$Final_Confidence_Grading, c(1L, 1L))
})
