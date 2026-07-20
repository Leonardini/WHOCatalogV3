make_comment_aux <- function(loF = NULL, category = NULL, single = NULL) {
  list(
    commentLoF         = loF      %||% tibble(drug = character(), gene = character(),
                                              comment = character()),
    commentCategoryTab = category %||% tibble(drug = character(), gene = character(),
                                              Final = integer(), comment = character(),
                                              exceptEffect = character()),
    commentSingleTab   = single   %||% tibble(drug = character(), gene = character(),
                                              mutation = character(), comment = character())
  )
}

test_that("applyComments backfills gene and mutation for variants added by expert rules", {
  inputTab <- tibble(drug = "Cycloserine", variant = "cycA_p.Gly122Ser",
                     gene = NA_character_, mutation = NA_character_,
                     effect_ALL = NA_character_, Final = 2L)
  auxData <- make_comment_aux(single = tibble(drug = "Cycloserine", gene = "cycA",
                                              mutation = "p.Gly122Ser",
                                              comment = "M. bovis BCG note"))
  out <- applyComments(inputTab, auxData)
  expect_equal(out$gene, "cycA")
  expect_equal(out$mutation, "p.Gly122Ser")
  expect_equal(out$comment, "M. bovis BCG note")
})

test_that("applyComments leaves an already parsed gene and mutation untouched", {
  inputTab <- tibble(drug = "Isoniazid", variant = "katG_p.Ser315Thr",
                     gene = "katG", mutation = "p.Ser315Thr",
                     effect_ALL = "missense_variant", Final = 1L)
  out <- applyComments(inputTab, make_comment_aux())
  expect_equal(out$gene, "katG")
  expect_equal(out$mutation, "p.Ser315Thr")
})

test_that("applyComments attaches the LoF comment to the pooled LoF row", {
  inputTab <- tibble(drug = "Bedaquiline", variant = "mmpL5_LoF", gene = "mmpL5",
                     mutation = "LoF", effect_ALL = "LoF", Final = 3L)
  auxData <- make_comment_aux(loF = tibble(drug = "Bedaquiline", gene = "mmpL5",
                                           comment = "Abrogates linked Rv0678 mutations"))
  out <- applyComments(inputTab, auxData)
  expect_equal(out$comment, "Abrogates linked Rv0678 mutations")
})

test_that("applyComments attaches the LoF comment to a constituent LoF effect", {
  inputTab <- tibble(drug = "Bedaquiline", variant = "mmpL5_p.Trp10*", gene = "mmpL5",
                     mutation = "p.Trp10*", effect_ALL = "stop_gained", Final = 3L)
  auxData <- make_comment_aux(loF = tibble(drug = "Bedaquiline", gene = "mmpL5",
                                           comment = "Abrogates linked Rv0678 mutations"))
  out <- applyComments(inputTab, auxData)
  expect_equal(out$comment, "Abrogates linked Rv0678 mutations")
})

test_that("applyComments withholds the LoF comment from a non-LoF effect", {
  inputTab <- tibble(drug = "Bedaquiline", variant = "mmpL5_p.Ala10Val", gene = "mmpL5",
                     mutation = "p.Ala10Val", effect_ALL = "missense_variant", Final = 3L)
  auxData <- make_comment_aux(loF = tibble(drug = "Bedaquiline", gene = "mmpL5",
                                           comment = "Abrogates linked Rv0678 mutations"))
  out <- applyComments(inputTab, auxData)
  expect_true(is.na(out$comment))
})

test_that("applyComments matches a category comment on the corrected grade", {
  auxData <- make_comment_aux(category = tibble(drug = "Moxifloxacin", gene = "gyrB",
                                                Final = 2L, comment = "Low-level resistance",
                                                exceptEffect = NA_character_))
  graded6 <- tibble(drug = "Moxifloxacin", variant = "gyrB_p.Asn499Thr", gene = "gyrB",
                    mutation = "p.Asn499Thr", effect_ALL = "missense_variant", Final = 6L)
  expect_true(is.na(applyComments(graded6, auxData)$comment))
  corrected <- graded6 %>% mutate(Final = 2L)
  expect_equal(applyComments(corrected, auxData)$comment, "Low-level resistance")
})
