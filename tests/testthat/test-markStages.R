make_stage_row <- function(tier, effect) {
  tibble(tier = tier, effect = effect)
}

test_that("markStages assigns stage 1 for tier 1 non-silent effect", {
  out <- markStages(make_stage_row(1L, "missense_variant"))
  expect_equal(out$stage, 1L)
})

test_that("markStages assigns stage 2 for tier 2 non-silent effect", {
  out <- markStages(make_stage_row(2L, "missense_variant"))
  expect_equal(out$stage, 2L)
})

test_that("markStages assigns stage 3 for tier 1 silent effect", {
  out <- markStages(make_stage_row(1L, "synonymous_variant"))
  expect_equal(out$stage, 3L)
})

test_that("markStages assigns stage 4 for tier 2 silent effect", {
  out <- markStages(make_stage_row(2L, "synonymous_variant"))
  expect_equal(out$stage, 4L)
})

test_that("markStages assigns stage 4 as default when tier is not 1 or 2", {
  out <- markStages(make_stage_row(3L, "missense_variant"))
  expect_equal(out$stage, 4L)
})

test_that("markStages assigns stage 3 for all three SILENT_EFFECTS at tier 1", {
  out <- markStages(tibble(tier = 1L, effect = SILENT_EFFECTS))
  expect_true(all(out$stage == 3L))
})
