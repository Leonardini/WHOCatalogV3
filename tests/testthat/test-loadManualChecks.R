testthat::local_edition(3)

nbsp_bytes <- function(...) {
  segs <- list(...)
  parts <- lapply(segs, function(x) if (identical(x, "NBSP")) as.raw(0xa0) else charToRaw(x))
  do.call(c, parts)
}

write_manual_check_fixture <- function(rowBytesList) {
  dir <- withr::local_tempdir(.local_envir = parent.frame())
  header <- charToRaw("drug,Gene,variant,Supplementary grading considerations,FINAL CONFIDENCE GRADING,Comment")
  body <- do.call(c, lapply(rowBytesList, function(rb) c(rb, as.raw(0x0a))))
  con <- file(file.path(dir, "manual_check.csv"), open = "wb")
  writeBin(c(header, as.raw(0x0a), body), con)
  close(con)
  dir
}

test_that("loadManualChecks removes internal non-breaking spaces so the variant can be matched", {
  cycloserineRow <- nbsp_bytes(
    "Cycloserine", "NBSP", ",cycA", "NBSP", ",cycA_", "NBSP", "p.Arg93Leu", "NBSP",
    ",Additional", "NBSP", "grading evidence: Literature evidence (PMID: 27064254; 28971867)",
    ",4) Not", "NBSP", "assoc", "NBSP", "w R - Interim,")
  dir <- write_manual_check_fixture(list(cycloserineRow))
  out <- loadManualChecks(dir)

  expect_equal(out$drug, "Cycloserine")
  expect_equal(out$variant, "cycA_p.Arg93Leu")
  expect_equal(out$Final_Confidence_Grading, "4) Not assoc w R - Interim")
  expect_equal(out$Supplementary_Grading_Considerations, "Evidence: Literature evidence (PMID: 27064254; 28971867)")
  expect_false(any(grepl("[[:space:]]", out$variant)))
})

test_that("loadManualChecks leaves variants without stray whitespace untouched", {
  cleanRow <- nbsp_bytes("Ethambutol,embB,embB_p.Gly406Asp,Additional grading evidence: WHO,2) Assoc w R - Interim,")
  dir <- write_manual_check_fixture(list(cleanRow))
  out <- loadManualChecks(dir)

  expect_equal(out$variant, "embB_p.Gly406Asp")
  expect_equal(out$Final_Confidence_Grading, "2) Assoc w R - Interim")
})
