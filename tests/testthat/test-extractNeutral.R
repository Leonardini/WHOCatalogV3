make_neutral_df <- function(set_col, values = c(TRUE, FALSE)) {
  tibble(
    drug    = c("Isoniazid", "Isoniazid"),
    variant = c("katG_p.Ser315Thr", "katG_p.Arg463Leu"),
    !!set_col := values
  )
}

test_that("extractNeutral returns only drug and variant columns for non-F sets", {
  out <- extractNeutral(make_neutral_df("setX"), "X")
  expect_equal(colnames(out), c("drug", "variant"))
})

test_that("extractNeutral excludes rows where the set column is FALSE", {
  out <- extractNeutral(make_neutral_df("setX"), "X")
  expect_equal(nrow(out), 1L)
  expect_equal(out$variant, "katG_p.Ser315Thr")
})

test_that("extractNeutral deduplicates rows with the same drug and variant", {
  df <- bind_rows(
    tibble(drug = "Isoniazid", variant = "katG_p.Ser315Thr", setX = TRUE),
    tibble(drug = "Isoniazid", variant = "katG_p.Ser315Thr", setX = TRUE)
  )
  out <- extractNeutral(df, "X")
  expect_equal(nrow(out), 1L)
})

test_that("extractNeutral for set F returns drug, variant, and all set membership columns", {
  df <- tibble(
    drug         = "Isoniazid",
    variant      = "katG_p.Ser315Thr",
    setF         = TRUE,
    setA         = TRUE,
    setB         = FALSE,
    setC         = TRUE,
    setD         = FALSE,
    setE         = FALSE,
    prev_version = FALSE,
    lit_mutation = FALSE
  )
  out <- extractNeutral(df, "F")
  expect_equal(
    colnames(out),
    c("drug", "variant", "setA", "setB", "setC", "setD", "setE", "prev_version", "lit_mutation")
  )
})

test_that("writeNeutralOutputs writes a CSV named neutral_mutations_WHO_{setName}.csv for each set", {
  td <- withr::local_tempdir()
  withr::local_dir(td)
  masterTab <- tibble(
    drug = "Isoniazid", variant = "katG_p.Ser315Thr",
    setA = TRUE, setB = TRUE, setC = TRUE, setD = TRUE, setE = TRUE, setF = TRUE,
    lit_mutation = FALSE, prev_version = FALSE
  )
  writeNeutralOutputs(masterTab)
  for (set in c("A", "B", "C", "D", "E", "F")) {
    expect_true(file.exists(paste0("neutral_mutations_WHO_", set, ".csv")))
  }
})
