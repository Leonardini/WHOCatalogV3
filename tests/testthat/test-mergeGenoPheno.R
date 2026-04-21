test_that("mergeGenoPheno returns a list split by phenotype group", {
  geno  <- tibble(drug = "DrugA", sample_id = 1L, variant = "V1", het = FALSE)
  pheno <- tibble(drug = "DrugA", sample_id = 1L, phenotype = "R", category_phenotype = "ALL")
  out   <- mergeGenoPheno(geno, pheno)
  expect_true("MAIN" %in% names(out))
  expect_equal(nrow(out$MAIN), 1L)
  expect_equal(out$MAIN$variant, "V1")
  expect_false(out$MAIN$het)
})

test_that("mergeGenoPheno sets het = TRUE for samples with no matching genotype", {
  geno  <- tibble(drug = "DrugA", sample_id = 1L, variant = "V1", het = FALSE)
  pheno <- tibble(drug = "DrugA", sample_id = c(1L, 2L), phenotype = c("R", "S"),
                  category_phenotype = "ALL")
  out   <- mergeGenoPheno(geno, pheno)$MAIN
  no_geno <- out[is.na(out$variant), ]
  expect_equal(nrow(no_geno), 1L)
  expect_true(no_geno$het)
})

test_that("mergeGenoPheno excludes phenotypes with unrecognised category", {
  geno  <- tibble(drug = "DrugA", sample_id = 1L, variant = "V1", het = FALSE)
  pheno <- tibble(drug = "DrugA", sample_id = 1L, phenotype = "R", category_phenotype = "OTHER")
  out   <- mergeGenoPheno(geno, pheno)
  expect_equal(length(out), 0L)
})

test_that("mergeGenoPheno includes both ALL and WHO phenotypes in the MAIN group", {
  geno  <- tibble(drug = "DrugA", sample_id = 1L, variant = "V1", het = FALSE)
  pheno <- tibble(drug = "DrugA", sample_id = c(1L, 1L), phenotype = c("R", "R"),
                  category_phenotype = c("ALL", "WHO"))
  out   <- mergeGenoPheno(geno, pheno)
  expect_true("MAIN" %in% names(out))
  expect_equal(nrow(out$MAIN), 2L)
})
