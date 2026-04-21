make_final_tab <- function() {
  tibble(
    selected  = c(TRUE,  FALSE, TRUE,  FALSE),
    phenotype = c("R",   "S",   "S",   "R"),
    drug      = "Isoniazid",
    N         = c(10L,   8L,    2L,    3L)
  )
}

make_final_tabs_list <- function() {
  list(finalTabs = list(
    LineageAny_RIFAll_Regular_1 = make_final_tab(),
    LineageAny_RIFAll_Relaxed_1 = make_final_tab()
  ))
}

make_epi_tabs_list <- function() {
  list(epiTabs = list(
    AMI = tibble(
      Group     = c("A",  "A",  "B",  "B"),
      variant   = "eis_p.Gly93Asp",
      phenotype = c("R",  "S",  "R",  "S"),
      drug      = "Amikacin",
      N         = c(5L,   2L,   8L,   15L)
    )
  ))
}

test_that("postprocessTabs returns fullTab with TP/TN/FP/FN and PPV/Sens/Spec columns", {
  result <- postprocessTabs(make_final_tabs_list(), version = CURR_VERSION,
                             relaxed = FALSE, sameRIF = TRUE, minQ = NA)
  expect_true("fullTab" %in% names(result))
  expect_true(all(c("TP", "TN", "FP", "FN", "PPV", "Sens", "Spec") %in% colnames(result$fullTab)))
})

test_that("postprocessTabs computes correct TP/TN/FP/FN counts", {
  result <- postprocessTabs(make_final_tabs_list(), version = CURR_VERSION,
                             relaxed = FALSE, sameRIF = TRUE, minQ = NA)
  row <- result$fullTab[result$fullTab$drug == "Isoniazid", ][1L, ]
  expect_equal(row$TP, 10L)
  expect_equal(row$TN, 8L)
  expect_equal(row$FP, 2L)
  expect_equal(row$FN, 3L)
})

test_that("postprocessTabs keeps all groups when minQ is NA", {
  result <- postprocessTabs(make_final_tabs_list(), version = CURR_VERSION,
                             relaxed = FALSE, sameRIF = TRUE, minQ = NA)
  expect_true(any(grepl("Regular", result$fullTab$group)))
  expect_true(any(grepl("Relaxed", result$fullTab$group)))
})

test_that("postprocessTabs filters to Regular groups when minQ is set and relaxed = FALSE", {
  result <- postprocessTabs(make_final_tabs_list(), version = CURR_VERSION,
                             relaxed = FALSE, sameRIF = TRUE, minQ = 1000L)
  expect_true(all(grepl("Regular", result$fullTab$group)))
})

test_that("postprocessTabs filters to Relaxed groups when minQ is set and relaxed = TRUE", {
  result <- postprocessTabs(make_final_tabs_list(), version = CURR_VERSION,
                             relaxed = TRUE, sameRIF = TRUE, minQ = 1000L)
  expect_true(all(grepl("Relaxed", result$fullTab$group)))
})

test_that("postprocessTabs returns no epiTab when only finalTabs is provided", {
  result <- postprocessTabs(make_final_tabs_list(), version = CURR_VERSION,
                             relaxed = FALSE, sameRIF = TRUE, minQ = NA)
  expect_false("epiTab" %in% names(result))
})

test_that("postprocessTabs returns epiTab with PPV_LOF and PPV_nomut columns", {
  result <- postprocessTabs(make_epi_tabs_list(), version = CURR_VERSION,
                             relaxed = FALSE, sameRIF = TRUE, minQ = NA)
  expect_true("epiTab" %in% names(result))
  expect_true(all(c("PPV_LOF", "PPV_nomut", "criterion") %in% colnames(result$epiTab)))
})

test_that("postprocessTabs epiTab computes correct LOF and nomut counts", {
  result <- postprocessTabs(make_epi_tabs_list(), version = CURR_VERSION,
                             relaxed = FALSE, sameRIF = TRUE, minQ = NA)
  row <- result$epiTab[result$epiTab$criterion == "eis_p.Gly93Asp", ][1L, ]
  expect_equal(row$LOF_R,   5L)
  expect_equal(row$LOF_S,   2L)
  expect_equal(row$nomut_R, 8L)
  expect_equal(row$nomut_S, 15L)
})

test_that("postprocessTabs fullTab has lb/ub CI columns for each metric", {
  result <- postprocessTabs(make_final_tabs_list(), version = CURR_VERSION,
                             relaxed = FALSE, sameRIF = TRUE, minQ = NA)
  expect_true(all(c("PPV_lb", "PPV_ub", "Sens_lb", "Sens_ub", "Spec_lb", "Spec_ub") %in%
                    colnames(result$fullTab)))
})

test_that("postprocessTabs fullTab group column matches the finalTabs list entry name", {
  tabs <- list(finalTabs = list(MyGroupName = make_final_tab()))
  result <- postprocessTabs(tabs, version = CURR_VERSION,
                             relaxed = FALSE, sameRIF = TRUE, minQ = NA)
  expect_equal(unique(result$fullTab$group), "MyGroupName")
})

test_that("postprocessTabs handles all-R samples with TN = 0 and FP = 0", {
  tab <- tibble(
    selected  = c(TRUE,  FALSE),
    phenotype = c("R",   "R"),
    drug      = "Isoniazid",
    N         = c(10L,   3L)
  )
  result <- postprocessTabs(list(finalTabs = list(test = tab)), version = CURR_VERSION,
                             relaxed = FALSE, sameRIF = TRUE, minQ = NA)
  row <- result$fullTab[1L, ]
  expect_equal(row$TN, 0L)
  expect_equal(row$FP, 0L)
})

test_that("postprocessTabs returns empty list for empty input", {
  result <- postprocessTabs(list(), version = CURR_VERSION,
                             relaxed = FALSE, sameRIF = TRUE, minQ = NA)
  expect_equal(length(result), 0L)
})
