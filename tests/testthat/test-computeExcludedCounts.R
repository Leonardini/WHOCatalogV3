make_excluded_curset <- function() {
  tibble(
    sample_id          = c("s1", "s2"),
    drug               = "Isoniazid",
    variant            = "katG_p.Ser315Thr",
    gene               = "katG",
    phenotype          = c("S", "R"),
    het                = FALSE,
    effect             = "missense_variant",
    tier               = 1L,
    neutral            = FALSE,
    category_phenotype = "ALL",
    RDen               = 2L,
    SDen               = 1L,
    stage              = 1L,
    lit_mutation       = FALSE,
    prev_version       = FALSE,
    setA = FALSE, setB = FALSE, setC = FALSE, setD = FALSE, setE = FALSE,
    pos1 = 315L,
    RRDR_NON_SILENT    = FALSE
  )
}

make_excluded_samplesToExclude <- function() {
  tibble(
    sample_id          = "s1",
    drug               = "Isoniazid",
    variant            = "katG_p.Ser315Thr",
    gene               = "katG",
    phenotype          = "S",
    het                = FALSE,
    effect             = "missense_variant",
    tier               = 1L,
    neutral            = FALSE,
    category_phenotype = "ALL"
  )
}

test_that("computeExcludedCounts returns a named list with extraCounts and extraSolos", {
  result <- computeExcludedCounts(make_excluded_curset(), make_excluded_samplesToExclude(), "MAIN")
  expect_true(is.list(result))
  expect_true(all(c("extraCounts", "extraSolos") %in% names(result)))
})

test_that("computeExcludedCounts extraCounts has expected columns", {
  result <- computeExcludedCounts(make_excluded_curset(), make_excluded_samplesToExclude(), "MAIN")
  expect_true(all(c("drug", "variant", "present_R", "present_S") %in% colnames(result$extraCounts)))
})

test_that("computeExcludedCounts extraSolos has expected columns", {
  result <- computeExcludedCounts(make_excluded_curset(), make_excluded_samplesToExclude(), "MAIN")
  expect_true(all(c("drug", "variant", "Rcnt", "Scnt") %in% colnames(result$extraSolos)))
})

test_that("computeExcludedCounts counts S phenotype in present_S for BAD_VAR_DRUG_PAIRS match", {
  result <- computeExcludedCounts(make_excluded_curset(), make_excluded_samplesToExclude(), "MAIN")
  expect_equal(result$extraCounts$present_S, 1L)
  expect_equal(result$extraCounts$present_R, 0L)
})

test_that("computeExcludedCounts returns empty extraCounts when no BAD_VAR_DRUG_PAIRS match", {
  curset <- make_excluded_curset() %>% mutate(drug = "Pyrazinamide", variant = "pncA_p.Trp119Stop", gene = "pncA")
  exclude <- make_excluded_samplesToExclude() %>% mutate(drug = "Pyrazinamide", variant = "pncA_p.Trp119Stop", gene = "pncA")
  result <- computeExcludedCounts(curset, exclude, "MAIN")
  expect_equal(nrow(result$extraCounts), 0L)
})

test_that("computeExcludedCounts with datasetName = 'MAIN' includes all samplesToExclude entries", {
  exclude <- bind_rows(
    make_excluded_samplesToExclude() %>% mutate(category_phenotype = "WHO"),
    make_excluded_samplesToExclude() %>% mutate(sample_id = "s3", category_phenotype = "ALL")
  )
  result <- computeExcludedCounts(make_excluded_curset(), exclude, "MAIN")
  expect_equal(result$extraCounts$present_S, 2L)
})

test_that("computeExcludedCounts returns empty outputs when nothing is excluded", {
  exclude <- make_excluded_samplesToExclude() %>% mutate(sample_id = "s99", category_phenotype = "ALL")
  result  <- computeExcludedCounts(make_excluded_curset(), exclude, "WHO")
  expect_equal(nrow(result$extraCounts), 0L)
  expect_equal(nrow(result$extraSolos),  0L)
})

test_that("computeExcludedCounts with datasetName = 'WHO' filters samplesToExclude by category_phenotype", {
  exclude <- bind_rows(
    make_excluded_samplesToExclude() %>% mutate(category_phenotype = "WHO"),
    make_excluded_samplesToExclude() %>% mutate(sample_id = "s3", category_phenotype = "ALL")
  )
  result <- computeExcludedCounts(make_excluded_curset(), exclude, "WHO")
  expect_equal(result$extraCounts$present_S, 1L)
})

