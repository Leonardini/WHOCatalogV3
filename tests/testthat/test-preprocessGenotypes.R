make_raw_genotypes <- function(af = 0.9, quality = 2000L,
                               mutation = "p.Ser315Thr", effect = "missense_variant") {
  tibble(
    drug             = "Isoniazid",
    sample_id        = "s1",
    resolved_symbol  = "katG",
    variant_category = mutation,
    predicted_effect = effect,
    neutral          = FALSE,
    tier             = 1L,
    `max(af)`        = af,
    `max(quality)`   = quality
  )
}

test_that("preprocessGenotypes renames resolved_symbol, variant_category, predicted_effect", {
  result <- preprocessGenotypes(make_raw_genotypes())
  expect_true(all(c("gene", "mutation", "effect") %in% colnames(result)))
  expect_false(any(c("resolved_symbol", "variant_category", "predicted_effect") %in% colnames(result)))
})

test_that("preprocessGenotypes drops the neutral column", {
  result <- preprocessGenotypes(make_raw_genotypes())
  expect_false("neutral" %in% colnames(result))
})

test_that("preprocessGenotypes creates variant as gene_mutation", {
  result <- preprocessGenotypes(make_raw_genotypes())
  expect_equal(result$variant, "katG_p.Ser315Thr")
})

test_that("preprocessGenotypes sets variant to 'missing' when mutation is NA", {
  result <- preprocessGenotypes(make_raw_genotypes(mutation = NA_character_, effect = NA_character_))
  expect_equal(result$variant, "missing")
})

test_that("preprocessGenotypes adds het = FALSE for high MAF and quality", {
  result <- preprocessGenotypes(make_raw_genotypes(af = 0.9, quality = 2000L))
  expect_false(result$het)
})

test_that("preprocessGenotypes marks het = TRUE for low MAF", {
  result <- preprocessGenotypes(make_raw_genotypes(af = 0.1))
  expect_true(result$het)
})

test_that("preprocessGenotypes adds pos1 column via extractPositions", {
  result <- preprocessGenotypes(make_raw_genotypes())
  expect_true("pos1" %in% colnames(result))
  expect_equal(result$pos1, 315L)
})
