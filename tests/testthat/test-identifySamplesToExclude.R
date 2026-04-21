make_main_dataset <- function() {
  tibble(
    sample_id          = c("s1", "s2", "s3", "s4"),
    drug               = c("Isoniazid", "Isoniazid", "Isoniazid", "Pyrazinamide"),
    variant            = c("katG_p.Ser315Thr", "katG_p.Ser315Thr", "katG_p.Arg463Leu", "pncA_p.Trp119Stop"),
    gene               = c("katG", "katG", "katG", "pncA"),
    phenotype          = c("S",    "R",    "S",    "S"),
    het                = FALSE,
    effect             = "missense_variant",
    tier               = 1L,
    category_phenotype = "ALL"
  )
}

test_that("identifySamplesToExclude returns only S-phenotype samples matching BAD_VAR_DRUG_PAIRS", {
  result <- identifySamplesToExclude(make_main_dataset())
  expect_equal(nrow(result), 1L)
  expect_equal(result$sample_id, "s1")
})

test_that("identifySamplesToExclude excludes R-phenotype samples even if variant matches", {
  result <- identifySamplesToExclude(make_main_dataset())
  expect_false("s2" %in% result$sample_id)
})

test_that("identifySamplesToExclude excludes samples with variant not in BAD_VAR_DRUG_PAIRS", {
  result <- identifySamplesToExclude(make_main_dataset())
  expect_false("s3" %in% result$sample_id)
  expect_false("s4" %in% result$sample_id)
})

test_that("identifySamplesToExclude sets neutral = FALSE for all returned rows", {
  result <- identifySamplesToExclude(make_main_dataset())
  expect_true(all(result$neutral == FALSE))
})

test_that("identifySamplesToExclude returns expected columns", {
  result <- identifySamplesToExclude(make_main_dataset())
  expect_true(all(c("sample_id", "drug", "variant", "gene", "phenotype", "het",
                    "effect", "tier", "neutral", "category_phenotype") %in% colnames(result)))
})
