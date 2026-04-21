test_that("describePMIDs formats a single PMID", {
  expect_equal(describePMIDs(12345678), make_Evidence("Literature evidence (PMID 12345678)"))
})

test_that("describePMIDs formats multiple PMIDs with semicolons", {
  expect_equal(describePMIDs(c(12345678, 87654321)), make_Evidence("Literature evidence (PMID 12345678; 87654321)"))
})
