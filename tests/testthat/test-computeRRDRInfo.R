test_that("returns data unchanged when RRDR_NON_SILENT column already exists", {
  df <- tibble(mutation = "p.Ser450Leu", gene = RRDR_GENE, effect = "missense_variant",
               RRDR_NON_SILENT = FALSE)
  out <- computeRRDRInfo(df)
  expect_equal(out$RRDR_NON_SILENT, FALSE)
  expect_false("pos1" %in% colnames(out))
})

test_that("sets RRDR_NON_SILENT = TRUE for RRDR gene, position in RRDR, non-silent effect", {
  df <- tibble(mutation = "p.Ser450Leu", gene = RRDR_GENE, effect = "missense_variant")
  out <- computeRRDRInfo(df)
  expect_true(out$RRDR_NON_SILENT)
})

test_that("sets RRDR_NON_SILENT = FALSE for a non-RRDR gene", {
  df <- tibble(mutation = "p.Ser450Leu", gene = "katG", effect = "missense_variant")
  out <- computeRRDRInfo(df)
  expect_false(out$RRDR_NON_SILENT)
})

test_that("sets RRDR_NON_SILENT = FALSE when position is outside the RRDR interval", {
  df <- tibble(mutation = "p.Ser531Leu", gene = RRDR_GENE, effect = "missense_variant")
  out <- computeRRDRInfo(df)
  expect_false(out$RRDR_NON_SILENT)
})

test_that("sets RRDR_NON_SILENT = FALSE for a silent effect in the RRDR gene and position", {
  df <- tibble(mutation = "p.Ser450Ser", gene = RRDR_GENE, effect = "synonymous_variant")
  out <- computeRRDRInfo(df)
  expect_false(out$RRDR_NON_SILENT)
})

test_that("uses existing pos1/pos2 columns without re-extracting positions from mutation", {
  df <- tibble(mutation = "not_a_real_mutation", gene = RRDR_GENE, effect = "missense_variant",
               pos1 = 450L, pos2 = NA_integer_)
  out <- computeRRDRInfo(df)
  expect_true(out$RRDR_NON_SILENT)
})

test_that("adds pos1 and pos2 columns when they are absent", {
  df <- tibble(mutation = "p.Ser450Leu", gene = RRDR_GENE, effect = "missense_variant")
  out <- computeRRDRInfo(df)
  expect_true("pos1" %in% colnames(out))
  expect_true("pos2" %in% colnames(out))
})
