make_annotate_dataset <- function() {
  df <- tibble(
    sample_id = c("s1", "s2", "s3"),
    drug      = "Pyrazinamide",
    variant   = "pncA_p.Trp119Stop",
    gene      = "pncA",
    phenotype = c("R", "S", "R"),
    het       = FALSE,
    category_phenotype = "ALL"
  )
  list(MAIN = df, WHO = df %>% mutate(category_phenotype = "WHO"))
}

make_exclude_one <- function() {
  tibble(sample_id = "s3", drug = "Pyrazinamide", variant = "pncA_p.Trp119Stop",
         phenotype = "R", het = FALSE, category_phenotype = "ALL")
}

make_empty_neutrals <- function() {
  tibble(drug = character(), variant = character(), neutral = logical(),
         setA = logical(), setB = logical(), setC = logical(),
         setD = logical(), setE = logical(),
         lit_mutation = logical(), prev_version = logical())
}

test_that("annotateDatasets adds RDen and SDen columns", {
  result <- annotateDatasets(make_annotate_dataset(), make_exclude_one(), make_empty_neutrals())
  expect_true(all(c("RDen", "SDen") %in% colnames(result$MAIN)))
})

test_that("annotateDatasets computes RDen and SDen excluding the excluded sample", {
  result <- annotateDatasets(make_annotate_dataset(), make_exclude_one(), make_empty_neutrals())
  row <- result$MAIN[1L, ]
  expect_equal(row$RDen, 1L)
  expect_equal(row$SDen, 1L)
})

test_that("annotateDatasets adds neutral set columns defaulting to FALSE", {
  result <- annotateDatasets(make_annotate_dataset(), make_exclude_one(), make_empty_neutrals())
  expect_true(all(c("setA", "setB", "neutral", "lit_mutation", "prev_version") %in% colnames(result$MAIN)))
  expect_true(all(result$MAIN$neutral == FALSE))
  expect_true(all(result$MAIN$setA    == FALSE))
})

test_that("annotateDatasets marks neutral = TRUE for variants in allNeutrals", {
  neutrals <- tibble(drug = "Pyrazinamide", variant = "pncA_p.Trp119Stop", neutral = TRUE,
                     setA = TRUE, setB = FALSE, setC = FALSE, setD = FALSE, setE = FALSE,
                     lit_mutation = FALSE, prev_version = FALSE)
  result <- annotateDatasets(make_annotate_dataset(), make_exclude_one(), neutrals)
  expect_true(all(result$MAIN$neutral == TRUE))
  expect_true(all(result$MAIN$setA    == TRUE))
})

test_that("annotateDatasets processes both MAIN and WHO entries", {
  result <- annotateDatasets(make_annotate_dataset(), make_exclude_one(), make_empty_neutrals())
  expect_true("RDen" %in% colnames(result$WHO))
})

test_that("annotateDatasets includes excluded bad-pair counts in denominators", {
  ds <- tibble(
    sample_id = c("s1", "s2", "s3"),
    drug      = "Isoniazid",
    variant   = "katG_p.Ser315Thr",
    gene      = "katG",
    phenotype = c("R", "S", "R"),
    het       = FALSE,
    category_phenotype = "ALL"
  )
  excluded <- tibble(sample_id = "s3", drug = "Isoniazid", variant = "katG_p.Ser315Thr",
                     phenotype = "R", het = FALSE, category_phenotype = "ALL")
  result <- annotateDatasets(list(MAIN = ds), excluded, make_empty_neutrals())
  katg_row <- result$MAIN[result$MAIN$variant == "katG_p.Ser315Thr", ][1L, ]
  expect_equal(katg_row$RDen, 2L)
  expect_equal(katg_row$SDen, 1L)
})
