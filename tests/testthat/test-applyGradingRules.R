make_grading_input <- function(who_SOLO_SorR = 0L, all_SOLO_SorR = 0L,
                                all_PPVc_SOLO_lb = NA_real_, all_OR_SOLO = NA_real_,
                                all_OR_SOLO_pval_FDR_sig = FALSE,
                                effect = "missense_variant") {
  base <- tibble(
    drug                 = "Isoniazid",
    variant              = "katG_p.Ser315Thr",
    gene                 = "katG",
    mutation             = "p.Ser315Thr",
    neutral              = NA,
    stage                = 2L,
    SOLO_SorR            = 0L,
    PPVc_SOLO_lb         = NA_real_,
    OR_SOLO              = NA_real_,
    OR_SOLO_pval_FDR_sig = FALSE,
    PPV_SOLO             = NA_real_,
    PPV_SOLO_ub          = NA_real_,
    SOLO_R               = 0L,
    PPVc_SOLO            = NA_real_,
    effect               = effect,
    pos1                 = 315L,
    pos2                 = NA_real_,
    setA                 = FALSE,
    setB                 = FALSE,
    setC                 = FALSE,
    setD                 = FALSE,
    setE                 = FALSE,
    prev_version         = FALSE
  )
  who_row <- base %>% mutate(datasets = "WHO", SOLO_SorR = who_SOLO_SorR)
  all_row <- base %>% mutate(
    datasets             = "ALL",
    SOLO_SorR            = all_SOLO_SorR,
    PPVc_SOLO_lb         = all_PPVc_SOLO_lb,
    OR_SOLO              = all_OR_SOLO,
    OR_SOLO_pval_FDR_sig = all_OR_SOLO_pval_FDR_sig
  )
  bind_rows(who_row, all_row)
}

make_aux_data <- function() {
  list(
    assayTab           = tibble(drug = character(), variant = character(), assay = logical()),
    allelicTab         = tibble(drug = character(), variant = character(), allelic = logical()),
    extraTab           = tibble(drug = character(), variant = character(),
                                final_grading_prev_version = character(),
                                Final_prev_version = integer()),
    commentCategoryTab = tibble(drug = character(), gene = character(),
                                Final = integer(), comment = character(),
                                exceptEffect = character()),
    commentLoF           = tibble(drug = character(), gene = character(), comment = character()),
    commentSingleTab     = tibble(drug = character(), gene = character(),
                                  mutation = character(), comment = character())
  )
}

test_that("applyGradingRules returns uncertain significance when no rules apply", {
  out <- applyGradingRules(make_grading_input(), make_aux_data())
  row <- out[out$drug == "Isoniazid", ]
  expect_equal(nrow(row), 1L)
  expect_equal(row$Final, GRADES[3L])
})

test_that("applyGradingRules applies the AllOnly rule when ALL has grade 1 but WHO has grade 3", {
  inputTab <- make_grading_input(
    all_SOLO_SorR            = 10L,
    all_PPVc_SOLO_lb         = 0.5,
    all_OR_SOLO              = 2,
    all_OR_SOLO_pval_FDR_sig = TRUE
  )
  out <- applyGradingRules(inputTab, make_aux_data())
  row <- out[out$drug == "Isoniazid", ]
  expect_equal(nrow(row), 1L)
  expect_equal(row$Final, GRADES[2L])
  expect_equal(row$`Additional grading criteria`, ALL_ONLY)
})

make_katG_comment_aux <- function() {
  aux <- make_aux_data()
  aux$commentCategoryTab <- tibble(
    drug         = "Isoniazid",
    gene         = "katG",
    Final        = 3L,
    comment      = "If genetically linked with an ahpC compensatory mutation, upgrade to Group 2",
    exceptEffect = "upstream_gene_variant"
  )
  aux
}

test_that("applyGradingRules suppresses a category comment for the excepted effect", {
  out <- applyGradingRules(make_grading_input(effect = "upstream_gene_variant"),
                           make_katG_comment_aux())
  row <- out[out$drug == "Isoniazid", ]
  expect_equal(nrow(row), 1L)
  expect_true(is.na(row$comment))
})

test_that("applyGradingRules keeps a category comment for effects that are not excepted", {
  out <- applyGradingRules(make_grading_input(), make_katG_comment_aux())
  row <- out[out$drug == "Isoniazid", ]
  expect_equal(nrow(row), 1L)
  expect_equal(row$comment,
               "If genetically linked with an ahpC compensatory mutation, upgrade to Group 2")
})

test_that("loadGradingAuxData parses the except clause out of a category mutation", {
  auxData <- loadGradingAuxData(system.file("extdata", package = "SOLOport"))
  katGRow <- auxData$commentCategoryTab %>%
    dplyr::filter(drug == "Isoniazid", gene == "katG", Final == 3L)
  expect_equal(nrow(katGRow), 1L)
  expect_equal(katGRow$exceptEffect, "upstream_gene_variant")
  expect_false("exceptEffect" %in% colnames(auxData$commentSingleTab))
  expect_false("exceptEffect" %in% colnames(auxData$commentLoF))
})
