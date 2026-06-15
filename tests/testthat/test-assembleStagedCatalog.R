make_stage_row <- function(drug, variant, stage, present, final) {
  tibble(
    drug          = drug,
    variant       = variant,
    stage         = stage,
    Final         = final,
    present_ALL   = present,
    present_R_ALL = present,
    present_S_ALL = 0,
    absent_R_ALL  = 0,
    absent_S_ALL  = 100 - present,
    present_WHO   = present,
    present_R_WHO = present,
    present_S_WHO = 0,
    absent_R_WHO  = 0,
    absent_S_WHO  = 80 - present,
    SOLO_R_ALL    = present
  )
}

make_stage_catalogs <- function() {
  stage1 <- bind_rows(
    make_stage_row("D", "v1", 1L, 10, 1),
    make_stage_row("D", "v2", 2L, 8,  3),
    make_stage_row("D", "v3", 3L, 5,  3)
  )
  stage2 <- bind_rows(
    make_stage_row("D", "v1", 1L, 10, 1),
    make_stage_row("D", "v2", 2L, 3,  1),
    make_stage_row("D", "v3", 3L, 4,  3)
  )
  stage3 <- bind_rows(
    make_stage_row("D", "v1", 1L, 10, 1),
    make_stage_row("D", "v2", 2L, 3,  1),
    make_stage_row("D", "v3", 3L, 1,  2),
    make_stage_row("D", "v4", NA_integer_, NA, NA)
  )
  list(stage1, stage2, stage3)
}

test_that("assembleStagedCatalog keeps one row per variant at its grading stage", {
  result <- assembleStagedCatalog(make_stage_catalogs())
  expect_equal(nrow(result), 4)
  expect_setequal(result$variant, c("v1", "v2", "v3", "v4"))
})

test_that("assembleStagedCatalog restores the true (stage 1) counts for late-stage variants", {
  result <- assembleStagedCatalog(make_stage_catalogs())
  expect_equal(result$present_ALL[result$variant == "v2"], 8)
  expect_equal(result$present_ALL[result$variant == "v3"], 5)
  expect_equal(result$present_ALL[result$variant == "v1"], 10)
})

test_that("assembleStagedCatalog restores every true-count column consistently", {
  result <- assembleStagedCatalog(make_stage_catalogs())
  v2 <- result[result$variant == "v2", ]
  expect_equal(v2$present_R_ALL, 8)
  expect_equal(v2$absent_S_ALL, 92)
  expect_equal(v2$present_WHO, 8)
  expect_equal(v2$absent_S_WHO, 72)
})

test_that("assembleStagedCatalog keeps grades and non-count columns from the grading stage", {
  result <- assembleStagedCatalog(make_stage_catalogs())
  expect_equal(result$Final[result$variant == "v2"], 1)
  expect_equal(result$Final[result$variant == "v3"], 2)
  expect_equal(result$SOLO_R_ALL[result$variant == "v3"], 1)
})

test_that("assembleStagedCatalog leaves ungraded (stage 4/NA) variants untouched", {
  result <- assembleStagedCatalog(make_stage_catalogs())
  expect_true(is.na(result$present_ALL[result$variant == "v4"]))
})
