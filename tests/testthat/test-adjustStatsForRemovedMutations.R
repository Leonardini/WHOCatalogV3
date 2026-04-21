make_cur_stats <- function() {
  tibble(
    drug        = "Isoniazid",
    variant     = "katG_p.Ser315Thr",
    present_R   = 3L,
    present_S   = 1L,
    absent_R    = 2L,
    absent_S    = 4L,
    SOLO_R      = 2L,
    SOLO_S      = 0L,
    SOLO_SorR   = 2L,
    stage       = 2L,
    tier        = 1L,
    neutral     = FALSE,
    datasets    = "ALL",
    effect      = "missense_variant",
    position    = "315",
    lit_mutation = FALSE,
    prev_version = FALSE,
    correctAll  = TRUE,
    correctSOLO = TRUE,
    pos1        = 315L
  )
}

make_cur_dataset <- function() {
  tibble(
    sample_id = c("s1", "s2", "s3", "s4"),
    drug      = "Isoniazid",
    variant   = "katG_p.Ser315Thr",
    phenotype = c("R", "R", "S", "R"),
    stage     = c(2L, 2L, 2L, 1L),
    het       = FALSE
  )
}

make_empty_solos <- function() {
  tibble(sample_id = character(), drug = character(), variant = character())
}

test_that("adjustStatsForRemovedMutations returns unchanged stats when no mutations match", {
  stats   <- make_cur_stats()
  dataset <- make_cur_dataset()
  remove  <- tibble(drug = "Rifampicin", variant = "rpoB_p.Ser450Leu")
  result  <- adjustStatsForRemovedMutations(stats, dataset, remove, 2L, make_empty_solos())
  expect_equal(result$present_R, stats$present_R)
  expect_equal(result$present_S, stats$present_S)
  expect_equal(result$absent_R,  stats$absent_R)
  expect_equal(result$absent_S,  stats$absent_S)
})

test_that("adjustStatsForRemovedMutations decrements present_R for removed R isolates", {
  stats   <- make_cur_stats()
  dataset <- make_cur_dataset()
  remove  <- tibble(drug = "Isoniazid", variant = "katG_p.Ser315Thr")
  result  <- adjustStatsForRemovedMutations(stats, dataset, remove, 2L, make_empty_solos())
  ## s1, s2 are R at stage 2; s3 is S at stage 2; s4 is R but at stage 1 (excluded)
  expect_equal(result$present_R, 1L)   ## 3 - 2 = 1
  expect_equal(result$present_S, 0L)   ## 1 - 1 = 0
})

test_that("adjustStatsForRemovedMutations recalculates absent_R and absent_S", {
  stats   <- make_cur_stats()
  dataset <- make_cur_dataset()
  remove  <- tibble(drug = "Isoniazid", variant = "katG_p.Ser315Thr")
  result  <- adjustStatsForRemovedMutations(stats, dataset, remove, 2L, make_empty_solos())
  ## RDen = present_R + absent_R = 3 + 2 = 5; new present_R = 1; absent_R = 5 - 1 = 4
  expect_equal(result$absent_R, 4L)
  ## SDen = present_S + absent_S = 1 + 4 = 5; new present_S = 0; absent_S = 5 - 0 = 5
  expect_equal(result$absent_S, 5L)
})

test_that("adjustStatsForRemovedMutations does not include RDen/SDen in output", {
  stats  <- make_cur_stats()
  result <- adjustStatsForRemovedMutations(stats, make_cur_dataset(),
                                           tibble(drug = "Rifampicin", variant = "rpoB_p.Ser450Leu"),
                                           2L, make_empty_solos())
  expect_false("RDen" %in% colnames(result))
  expect_false("SDen" %in% colnames(result))
})

test_that("adjustStatsForRemovedMutations subtracts SOLO counts for removed isolates", {
  stats   <- make_cur_stats()
  dataset <- make_cur_dataset()
  remove  <- tibble(drug = "Isoniazid", variant = "katG_p.Ser315Thr")
  solos   <- tibble(sample_id = "s1", drug = "Isoniazid", variant = "katG_p.Ser315Thr")
  result  <- adjustStatsForRemovedMutations(stats, dataset, remove, 2L, solos)
  expect_equal(result$SOLO_R, 1L)   ## 2 - 1 = 1
  expect_equal(result$SOLO_S, 0L)
  expect_equal(result$SOLO_SorR, 1L)
})

test_that("adjustStatsForRemovedMutations ignores solos for isolates not in extraIsolates", {
  stats   <- make_cur_stats()
  dataset <- make_cur_dataset()
  remove  <- tibble(drug = "Rifampicin", variant = "rpoB_p.Ser450Leu")
  solos   <- tibble(sample_id = "s1", drug = "Isoniazid", variant = "katG_p.Ser315Thr")
  result  <- adjustStatsForRemovedMutations(stats, dataset, remove, 2L, solos)
  expect_equal(result$SOLO_R, stats$SOLO_R)
})

test_that("adjustStatsForRemovedMutations prints a message", {
  stats   <- make_cur_stats()
  dataset <- make_cur_dataset()
  remove  <- tibble(drug = "Isoniazid", variant = "katG_p.Ser315Thr")
  expect_output(
    adjustStatsForRemovedMutations(stats, dataset, remove, 2L, make_empty_solos()),
    "Removing"
  )
})
