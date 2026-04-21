make_curset <- function(het = FALSE, variant = "katG_p.Ser315Thr",
                         effect = "missense_variant",
                         lit_mutation = FALSE, prev_version = FALSE,
                         rcnt = 2L, scnt = 0L) {
  tibble(
    drug              = "Isoniazid",
    variant           = variant,
    sample_id         = c("s1", "s2", "s3"),
    phenotype         = c("R", "R", "S"),
    het               = het,
    category_phenotype = "ALL",
    RDen              = 3L,
    SDen              = 1L,
    stage             = 1L,
    tier              = 1L,
    neutral           = FALSE,
    effect            = effect,
    position          = "315",
    lit_mutation      = lit_mutation,
    prev_version      = prev_version,
    Rcnt              = rcnt,
    Scnt              = scnt,
    pos1              = 315L,
    setA              = FALSE
  )
}

make_empty_extra_counts <- function() {
  tibble(drug = character(), variant = character(), present_R = integer(), present_S = integer())
}

make_empty_extra_solos <- function() {
  tibble(drug = character(), variant = character(), Rcnt = integer(), Scnt = integer())
}

test_that("computeVariantStats returns expected output columns", {
  result <- computeVariantStats(make_curset(), make_empty_extra_counts(), make_empty_extra_solos())
  expected_cols <- c("drug", "variant", "stage", "tier", "neutral", "datasets",
                     "effect", "position", "lit_mutation", "prev_version",
                     "present", "present_R", "present_S",
                     "absent_R", "absent_S",
                     "SOLO_R", "SOLO_S", "SOLO_SorR",
                     "correctAll", "correctSOLO", "pos1", "setA")
  expect_true(all(expected_cols %in% colnames(result)))
})

test_that("computeVariantStats renames category_phenotype to datasets", {
  result <- computeVariantStats(make_curset(), make_empty_extra_counts(), make_empty_extra_solos())
  expect_true("datasets" %in% colnames(result))
  expect_false("category_phenotype" %in% colnames(result))
})

test_that("computeVariantStats filters out het rows", {
  curset <- bind_rows(
    make_curset(het = FALSE),
    make_curset(het = TRUE) %>% mutate(sample_id = c("s4", "s5", "s6"))
  )
  result <- computeVariantStats(curset, make_empty_extra_counts(), make_empty_extra_solos())
  expect_equal(result$present, 3L)
})

test_that("computeVariantStats filters out missing variants", {
  curset <- bind_rows(
    make_curset(),
    make_curset(variant = "missing") %>% mutate(sample_id = c("s4", "s5", "s6"))
  )
  result <- computeVariantStats(curset, make_empty_extra_counts(), make_empty_extra_solos())
  expect_equal(nrow(result), 1L)
  expect_equal(result$variant, "katG_p.Ser315Thr")
})

test_that("computeVariantStats computes present_R and present_S", {
  result <- computeVariantStats(make_curset(), make_empty_extra_counts(), make_empty_extra_solos())
  expect_equal(result$present_R, 2L)
  expect_equal(result$present_S, 1L)
})

test_that("computeVariantStats computes absent_R and absent_S from RDen/SDen", {
  result <- computeVariantStats(make_curset(), make_empty_extra_counts(), make_empty_extra_solos())
  expect_equal(result$absent_R, 1L)  ## RDen=3, present_R=2
  expect_equal(result$absent_S, 0L)  ## SDen=1, present_S=1
})

test_that("computeVariantStats adds curExtraCounts to present counts", {
  extra <- tibble(drug = "Isoniazid", variant = "katG_p.Ser315Thr",
                  present_R = 1L, present_S = 0L)
  result <- computeVariantStats(make_curset(), extra, make_empty_extra_solos())
  expect_equal(result$present_R, 3L)
  expect_equal(result$present_S, 1L)
})

test_that("computeVariantStats adds curExtraSolos to SOLO counts", {
  extra <- tibble(drug = "Isoniazid", variant = "katG_p.Ser315Thr",
                  Rcnt = 1L, Scnt = 0L)
  result <- computeVariantStats(make_curset(), make_empty_extra_counts(), extra)
  expect_equal(result$SOLO_R, 3L)
})

test_that("computeVariantStats sets correctAll = FALSE for silent effect", {
  result <- computeVariantStats(
    make_curset(effect = "synonymous_variant"),
    make_empty_extra_counts(), make_empty_extra_solos()
  )
  expect_false(result$correctAll)
})

test_that("computeVariantStats sets correctAll = FALSE when lit_mutation = TRUE", {
  result <- computeVariantStats(
    make_curset(lit_mutation = TRUE),
    make_empty_extra_counts(), make_empty_extra_solos()
  )
  expect_false(result$correctAll)
})

test_that("computeVariantStats sets correctSOLO = FALSE when SOLO_SorR is NA", {
  result <- computeVariantStats(
    make_curset(rcnt = NA_integer_, scnt = NA_integer_),
    make_empty_extra_counts(), make_empty_extra_solos()
  )
  expect_false(result$correctSOLO)
})

test_that("computeVariantStats sets correctSOLO = FALSE when SOLO_SorR is 0", {
  result <- computeVariantStats(
    make_curset(rcnt = 0L, scnt = 0L),
    make_empty_extra_counts(), make_empty_extra_solos()
  )
  expect_false(result$correctSOLO)
})

test_that("computeVariantStats computes SOLO_SorR as SOLO_R + SOLO_S", {
  result <- computeVariantStats(
    make_curset(rcnt = 2L, scnt = 1L),
    make_empty_extra_counts(), make_empty_extra_solos()
  )
  expect_equal(result$SOLO_SorR, 3L)
})
