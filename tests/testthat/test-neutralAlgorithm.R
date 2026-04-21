make_neutral_input <- function(n_samples = 5L) {
  tibble(
    drug           = "Isoniazid",
    variant        = "katG_p.Arg463Leu",
    gene           = "katG",
    mutation       = "p.Arg463Leu",
    sample_id      = seq_len(n_samples),
    het            = FALSE,
    phenotype      = "S",
    effect         = "missense_variant",
    tier           = 1L,
    `max(af)`      = 1.0,
    `max(quality)` = 2000L
  )
}

make_catalogue_file <- function(td) {
  path <- file.path(td, "catalogue.csv")
  write_csv(
    tibble(drug = "DummyDrug", variant = "DummyGene_DummyMut",
           Final_Confidence_Grading = GRADES[3L]),
    path
  )
  path
}

make_lit_tab <- function(drug = "DummyDrug", variant = "DummyGene_DummyMut") {
  tibble(drug = drug, variant = variant, literature = 1L)
}

test_that("neutralAlgorithm returns the expected set columns", {
  td <- withr::local_tempdir()
  out <- neutralAlgorithm(
    make_neutral_input(),
    litTab        = make_lit_tab(),
    catalogueFile = make_catalogue_file(td)
  )
  expect_true(all(
    c("setA", "setB", "setC", "setD", "setE", "setF", "prev_version", "lit_mutation") %in%
      colnames(out)
  ))
})

test_that("neutralAlgorithm puts a variant with PPV_ub below threshold into setA and setF", {
  td <- withr::local_tempdir()
  out <- neutralAlgorithm(
    make_neutral_input(40L),
    litTab        = make_lit_tab(),
    catalogueFile = make_catalogue_file(td)
  )
  row <- out[out$drug == "Isoniazid" & out$variant == "katG_p.Arg463Leu", ][1L, ]
  expect_true(row$setA)
  expect_true(row$setF)
})

test_that("neutralAlgorithm puts a literature variant into setC and setF but not setA", {
  td <- withr::local_tempdir()
  out <- neutralAlgorithm(
    make_neutral_input(5L),
    litTab        = make_lit_tab(drug = "Isoniazid", variant = "katG_p.Arg463Leu"),
    catalogueFile = make_catalogue_file(td)
  )
  row <- out[out$drug == "Isoniazid" & out$variant == "katG_p.Arg463Leu", ][1L, ]
  expect_false(row$setA)
  expect_true(row$setC)
  expect_true(row$setF)
})
