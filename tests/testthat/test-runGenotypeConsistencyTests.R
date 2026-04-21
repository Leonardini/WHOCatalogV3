make_geno <- function(gene = "katG", drug = "Isoniazid",
                      mutation = "p.Ser315Thr", effect = "missense_variant",
                      pos1 = 315L, pos2 = NA_real_, position = 315L,
                      variant = "katG_p.Ser315Thr", tier = 1L) {
  tibble(gene = gene, drug = drug, mutation = mutation, effect = effect,
         pos1 = pos1, pos2 = pos2, position = position,
         variant = variant, tier = tier)
}

test_that("runGenotypeConsistencyTests passes for valid data", {
  expect_no_error(runGenotypeConsistencyTests(make_geno()))
})

test_that("runGenotypeConsistencyTests errors when pncA row has wrong drug", {
  testthat::local_edition(3)
  df <- make_geno(gene = "pncA", drug = "Isoniazid",
                  mutation = "p.Ala3Val", pos1 = 3L, position = 3L,
                  variant = "pncA_p.Ala3Val")
  expect_snapshot(runGenotypeConsistencyTests(df), error = TRUE)
})

test_that("runGenotypeConsistencyTests errors when rpoB row has wrong drug", {
  testthat::local_edition(3)
  df <- make_geno(gene = "rpoB", drug = "Isoniazid",
                  mutation = "p.Ser450Leu", pos1 = 450L, position = 450L,
                  variant = "rpoB_p.Ser450Leu")
  expect_snapshot(runGenotypeConsistencyTests(df), error = TRUE)
})

test_that("runGenotypeConsistencyTests errors when mutation is missing but effect is not", {
  testthat::local_edition(3)
  df <- make_geno(mutation = "missing", effect = "missense_variant",
                  pos1 = NA_integer_, position = NA_integer_,
                  variant = "katG_missing")
  expect_snapshot(runGenotypeConsistencyTests(df), error = TRUE)
})

test_that("runGenotypeConsistencyTests errors when effect is inconsistent for the same variant", {
  testthat::local_edition(3)
  df <- bind_rows(make_geno(), make_geno(effect = "stop_gained"))
  expect_snapshot(runGenotypeConsistencyTests(df), error = TRUE)
})

test_that("runGenotypeConsistencyTests errors when tier is inconsistent for the same drug-gene", {
  testthat::local_edition(3)
  df <- bind_rows(
    make_geno(),
    make_geno(mutation = "p.Arg463Leu", pos1 = 463L, position = 463L,
              variant = "katG_p.Arg463Leu", tier = 2L)
  )
  expect_snapshot(runGenotypeConsistencyTests(df), error = TRUE)
})
