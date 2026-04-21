make_pheno <- function(sample_id = 1L, category_phenotype = "ALL",
                       drug = "Isoniazid", phenotype = "R") {
  tibble(sample_id = sample_id, category_phenotype = category_phenotype,
         drug = drug, phenotype = phenotype)
}

test_that("runPhenotypeConsistencyTests passes for valid data", {
  expect_no_error(runPhenotypeConsistencyTests(make_pheno()))
})

test_that("runPhenotypeConsistencyTests errors when phenotype is inconsistent for same sample-category-drug", {
  testthat::local_edition(3)
  df <- bind_rows(make_pheno(phenotype = "R"), make_pheno(phenotype = "S"))
  expect_snapshot(runPhenotypeConsistencyTests(df), error = TRUE)
})

test_that("runPhenotypeConsistencyTests errors when same sample-drug appears in both WHO and ALL", {
  testthat::local_edition(3)
  df <- bind_rows(make_pheno(category_phenotype = "WHO"),
                  make_pheno(category_phenotype = "ALL"))
  expect_snapshot(runPhenotypeConsistencyTests(df), error = TRUE)
})
