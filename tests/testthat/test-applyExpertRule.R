make_expert_df <- function(anyRule = FALSE, ruleVal = TRUE, final = 3L, rule_final = 15L,
                            description = NA_character_) {
  tibble(
    anyRule               = anyRule,
    rule_test             = ruleVal,
    Final                 = final,
    Rule_Final            = rule_final,
    `Additional grading criteria` = description
  )
}

test_that("applyExpertRule sets Final, description and Rule_Final when rule applies", {
  out <- applyExpertRule(make_expert_df(), "rule_test",
                         description = "Upgraded", finalGrade = 2L, finalRule = 99L)
  expect_equal(out$Final, 2L)
  expect_equal(out$`Additional grading criteria`, "Upgraded")
  expect_equal(out$Rule_Final, 99L)
  expect_true(out$anyRule)
})

test_that("applyExpertRule does not apply when anyRule is already TRUE", {
  out <- applyExpertRule(make_expert_df(anyRule = TRUE), "rule_test",
                         description = "Upgraded", finalGrade = 2L, finalRule = 99L)
  expect_equal(out$Final, 3L)
  expect_true(is.na(out$`Additional grading criteria`))
  expect_equal(out$Rule_Final, 15L)
})

test_that("applyExpertRule with applyAlways = TRUE overrides anyRule", {
  out <- applyExpertRule(make_expert_df(anyRule = TRUE), "rule_test",
                         description = "Always", finalGrade = 1L, finalRule = 50L,
                         applyAlways = TRUE)
  expect_equal(out$Final, 1L)
  expect_equal(out$`Additional grading criteria`, "Always")
  expect_equal(out$Rule_Final, 50L)
})

test_that("applyExpertRule with finalGrade = NA leaves Final unchanged", {
  out <- applyExpertRule(make_expert_df(final = 2L), "rule_test",
                         description = "Note", finalGrade = NA, finalRule = 99L)
  expect_equal(out$Final, 2L)
  expect_equal(out$`Additional grading criteria`, "Note")
})

test_that("applyExpertRule with description = NA leaves Additional grading criteria unchanged", {
  out <- applyExpertRule(make_expert_df(description = "Previous"), "rule_test",
                         description = NA, finalGrade = 2L, finalRule = 99L)
  expect_equal(out$`Additional grading criteria`, "Previous")
})

test_that("applyExpertRule does not apply when ruleColumn is FALSE", {
  out <- applyExpertRule(make_expert_df(ruleVal = FALSE), "rule_test",
                         description = "Upgraded", finalGrade = 2L, finalRule = 99L)
  expect_equal(out$Final, 3L)
  expect_false(out$anyRule)
})
