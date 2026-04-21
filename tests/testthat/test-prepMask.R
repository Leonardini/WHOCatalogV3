make_mask_df <- function() {
  tibble(
    drug      = "DrugA",
    sample_id = 1:3,
    variant   = c("V1", "V2", "V3"),
    phenotype = c("R", "S", "S"),
    het       = FALSE,
    effect    = "missense_variant",
    tier      = 1L,
    neutral   = FALSE
  )
}

test_that("prepMask with SOnly = TRUE returns SOnly column and not het", {
  out <- prepMask(make_mask_df(), SOnly = TRUE)
  expect_true("SOnly" %in% colnames(out))
  expect_false("het" %in% colnames(out))
})

test_that("prepMask with SOnly = FALSE returns het column and not SOnly", {
  out <- prepMask(make_mask_df(), SOnly = FALSE)
  expect_true("het" %in% colnames(out))
  expect_false("SOnly" %in% colnames(out))
})

test_that("prepMask marks variants appearing only in S samples as SOnly = TRUE", {
  out <- prepMask(make_mask_df(), SOnly = TRUE)
  expect_false(out$SOnly[out$variant == "V1"])
  expect_true(out$SOnly[out$variant == "V2"])
  expect_true(out$SOnly[out$variant == "V3"])
})

test_that("prepMask removes tier 2 variants when Tier2 = TRUE", {
  df <- make_mask_df() %>% mutate(tier = c(1L, 2L, 1L))
  out <- prepMask(df, SOnly = FALSE, Tier2 = TRUE)
  expect_equal(nrow(out), 2L)
  expect_false("V2" %in% out$variant)
})

test_that("prepMask keeps tier 2 variants when Tier2 = FALSE", {
  df <- make_mask_df() %>% mutate(tier = c(1L, 2L, 1L))
  out <- prepMask(df, SOnly = FALSE, Tier2 = FALSE)
  expect_equal(nrow(out), 3L)
})

test_that("prepMask removes neutral variants when Neutral = TRUE", {
  df <- make_mask_df() %>% mutate(neutral = c(FALSE, TRUE, FALSE))
  out <- prepMask(df, SOnly = FALSE, Neutral = TRUE)
  expect_equal(nrow(out), 2L)
  expect_false("V2" %in% out$variant)
})

test_that("prepMask removes silent effects when Silent = TRUE", {
  df <- make_mask_df() %>% mutate(effect = c("missense_variant", "synonymous_variant", "missense_variant"))
  out <- prepMask(df, SOnly = FALSE, Silent = TRUE)
  expect_equal(nrow(out), 2L)
  expect_false("V2" %in% out$variant)
})

test_that("prepMask retains 'missing' effect even when Silent = TRUE", {
  df <- make_mask_df() %>% mutate(effect = "missing")
  out <- prepMask(df, SOnly = FALSE, Silent = TRUE)
  expect_equal(nrow(out), 3L)
})
