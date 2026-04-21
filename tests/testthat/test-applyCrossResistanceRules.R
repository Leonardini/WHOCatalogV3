make_cr_row <- function(drug, gene, effect = "missense_variant",
                         initial = 1L, final = 1L, any_rule = FALSE) {
  tibble(
    drug                        = drug,
    variant                     = paste0(gene, "_p.Ala90Val"),
    gene                        = gene,
    effect_ALL                  = effect,
    Initial                     = initial,
    Final                       = final,
    anyRule                     = any_rule,
    `Additional grading criteria` = NA_character_,
    Rule_Final                  = 15L
  )
}

test_that("applyCrossResistanceRules adds MXF row with grade 2 for LEV variant in gyrA", {
  out <- applyCrossResistanceRules(make_cr_row("Levofloxacin", "gyrA"),
                                   iteration = 1L, initRule = 0)
  mxf <- out[out$drug == "Moxifloxacin", ]
  expect_equal(nrow(mxf), 1L)
  expect_true(mxf$Final == 2)
  expect_equal(mxf$`Additional grading criteria`, CROSS_RES$description[1])
})

test_that("applyCrossResistanceRules does not add MXF for a variant in a non-FQ gene", {
  out <- applyCrossResistanceRules(make_cr_row("Levofloxacin", "katG"),
                                   iteration = 1L, initRule = 0)
  expect_equal(nrow(out[out$drug == "Moxifloxacin", ]), 0L)
})

test_that("applyCrossResistanceRules adds ETO row with grade 2 for INH variant in inhA", {
  out <- applyCrossResistanceRules(
    make_cr_row("Isoniazid", "inhA", effect = "upstream_gene_variant"),
    iteration = 1L, initRule = 0
  )
  eto <- out[out$drug == "Ethionamide", ]
  expect_equal(nrow(eto), 1L)
  expect_true(eto$Final == 2)
  expect_equal(eto$`Additional grading criteria`, CROSS_RES$description[2])
})

test_that("applyCrossResistanceRules adds CFZ row with grade 2 for BDQ variant in Rv0678", {
  out <- applyCrossResistanceRules(make_cr_row("Bedaquiline", "Rv0678"),
                                   iteration = 1L, initRule = 0)
  cfz <- out[out$drug == "Clofazimine", ]
  expect_equal(nrow(cfz), 1L)
  expect_true(cfz$Final == 2)
  expect_equal(cfz$`Additional grading criteria`, CROSS_RES$description[3])
})

test_that("applyCrossResistanceRules upgrades the mirror FQ row even when the source row has anyRule TRUE", {
  out <- applyCrossResistanceRules(
    make_cr_row("Levofloxacin", "gyrA", any_rule = TRUE),
    iteration = 1L, initRule = 0
  )
  mxf <- out[out$drug == "Moxifloxacin", ]
  expect_equal(nrow(mxf), 1L)
  expect_true(mxf$Final == 2)
})
