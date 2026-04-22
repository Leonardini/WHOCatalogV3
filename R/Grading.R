## TODO: fbiA_c.744C>A with PMD and DLM is changing from group 4 to group 2, confirm that this is correct!
## TODO: Take care of the inframe insertions and deletions and missense mutations (making sure they are listed for both) when updating!

GRADING_RULE_LEDGER = tribble(
  ~name,                  ~description,              ~grade,          ~applyAlways, ~nRules,
  "rule_literature",      LIT_EVIDENCE,              4L,              FALSE,        1L,
  "rule_AllOnly",         ALL_ONLY,                  2L,              FALSE,        1L,
  "rule_WHOBased",        NA_character_,             NA_integer_,     FALSE,        1L,
  "rule_FLAG",            MANUAL_CHECK,              NA_integer_,     FALSE,        1L,
  "rule_rpoB_borderline", NA_character_,             1L,              TRUE,         1L,
  "rule_silent",          SILENT,                    4L,              FALSE,        1L,
  "rule_RRDR",            RRDR,                      2L,              FALSE,        1L,
  "rule_LoF",             LoF_MUT,                   2L,              FALSE,        1L,
  "cross_res",            NA_character_,             2L,              FALSE,        NUM_CROSS_RES_RULES,
  "rule_Assay",           PROBLEMATIC_CRITERION,     2L,              FALSE,        1L,
  "rule_Allelic",         NA_character_,             2L,              FALSE,        1L,
  "rule_LoF_re",          LoF_MUT,                   2L,              FALSE,        1L,
  "rule_PZA_lit",         NA_character_,             2L,              FALSE,        1L,
  "rule_prev_guidance",   PREV_GUIDANCE,             FINAL_FLAG,      FALSE,        1L
) %>%
  mutate(.endRule = FIRST_RULE_NUM - 1L + cumsum(nRules))

ruleEndNum = function(name) GRADING_RULE_LEDGER$.endRule[GRADING_RULE_LEDGER$name == name]

applyNamedRule = function(inputTab, name, description = NULL, finalRule = NULL) {
  row  = GRADING_RULE_LEDGER[GRADING_RULE_LEDGER$name == name, ]
  desc = if (is.null(description)) row$description[[1]] else description
  fr   = if (is.null(finalRule))   row$.endRule[[1]]   else finalRule
  inputTab %>%
    applyExpertRule(name, description = desc, finalGrade = row$grade[[1]],
                    finalRule = fr, applyAlways = row$applyAlways[[1]])
}

#' Prepare the master variant table, carry out basic sanity checks, and compute genes and mutations
#' @noRd
prepareGradingInput = function(stage = NA, LoF = TRUE, outDir = ".") {
  Tab0 = read_csv(file.path(outDir, paste0(STATS_FILE_PREFIX, "WHO" , ifelse(LoF, WITHLOFS_SUFFIX, ""), "_Stage", stage, ".csv")), guess_max = LARGE_NUMBER, show_col_types = FALSE)
  Tab1 = read_csv(file.path(outDir, paste0(STATS_FILE_PREFIX, "MAIN", ifelse(LoF, WITHLOFS_SUFFIX, ""), "_Stage", stage, ".csv")), guess_max = LARGE_NUMBER, show_col_types = FALSE)
  stopifnot(all(Tab0$datasets == "WHO"))
  stopifnot(all(Tab1$datasets == "ALL"))
  stopifnot(nrow(Tab0 %>% anti_join(Tab1, by = c("variant", "drug"))) == 0)
  inputTab = bind_rows(Tab0, Tab1) %>%
    mutate(gene = str_split_fixed(variant, "_", n = 2)[,1], mutation = str_split_fixed(variant, "_", n = 2)[,2])
  inputTab
}

#' Load auxiliary grading data files (assay mutations, allelic exchanges, comments, and v2 grades)
#' @noRd
loadGradingAuxData = function(dir) {
  dir = str_remove(dir, "/$")
  assayTab = read_csv(paste0(dir, "/assay_mutations.csv"), guess_max = LARGE_NUMBER, show_col_types = FALSE) %>%
    mutate(drug = DRUG_LIST[match(drug, SHORT_NAMES)]) %>%
    mutate(assay = TRUE)
  allelicTab = read_csv(paste0(dir, "/allelic_exchanges.csv"), show_col_types = FALSE, guess_max = LARGE_NUMBER) %>%
    mutate(variant = str_replace(variant, "lof$", LOF_LABEL)) %>%
    mutate(allelic = TRUE)
  extraTab = read_csv(paste0(dir, "/v", PREV_VERSION, "_grades.csv"), guess_max = LARGE_NUMBER, show_col_types = FALSE) %>%
    mutate(Final_prev_version = match(final_grading_prev_version, GRADES))
  commentTab = read_csv(paste0(dir, "/grading_comments.csv"), show_col_types = FALSE, guess_max = LARGE_NUMBER) %>%
    mutate_all(~{trimws(., whitespace = "[\\h\\v]")}) %>%
    select(1:4) %>%
    set_colnames(c("drug", "gene", "mutation", "comment")) %>%
    mutate(temp = match(mutation, c("any AwR", "any AwRI", "any Uncertain", "any nAwRI", "any nAwR"))) %>%
    mutate(mutation = ifelse(!is.na(temp), temp, mutation)) %>%
    select(-temp)
  list(
    assayTab           = assayTab,
    allelicTab         = allelicTab,
    extraTab           = extraTab,
    commentCategoryTab = commentTab %>%
      dplyr::filter(nchar(mutation) == 1) %>%
      mutate(Final = as.integer(mutation)) %>%
      select(-mutation),
    commentLoF         = commentTab %>%
      dplyr::filter(str_starts(mutation, "any LoF")) %>%
      select(-mutation),
    commentSingleTab   = commentTab %>%
      dplyr::filter(nchar(mutation) > 1 & !str_starts(mutation, "any LoF"))
  )
}

#' Apply grading rules to all variants and write the graded catalogue
#' @param stage Integer stage index (1, 2, or 3); NA runs all stages together.
#' @param outDir Directory for reading intermediate stat files and writing grading output.
#' @inheritParams mainDriver
#' @return A data frame of graded variants with one row per variant-drug pair, including columns
#'   \code{variant}, \code{drug}, \code{Initial} and \code{Final} (grade labels from \code{GRADES}),
#'   \code{Additional grading criteria}, and rule-tracking columns. Also writes a dated CSV to \code{outDir}.
#' @export
gradeMutations = function(LoF = TRUE, NON_DATABASE_DIRECTORY = system.file("extdata", package = "SOLOport"), stage = NA, outDir = ".") {
  inputTab = prepareGradingInput(stage = stage, LoF = LoF, outDir = outDir)
  auxData = loadGradingAuxData(NON_DATABASE_DIRECTORY)
  inputTab = applyGradingRules(inputTab, auxData)
  writeGradingOutput(inputTab, stage = stage, LoF = LoF, outDir = outDir)
  inputTab
}

#' Apply the full set of grading rules to a variant table
#' @param inputTab A data frame of variants prepared by prepareGradingInput.
#' @param auxData A list of auxiliary grading data as returned by loadGradingAuxData.
#' @export
applyGradingRules = function(inputTab, auxData) {
  ## Initial confidence grading: upgrade variants matching the basic upgrade rule to grade 1; apply rules for upgrading and downgrading pncA variants
  ## NOTE: THE UPGRADE RULE HAS EXPANDED
  ## NOTE: The up and downgrade rules only apply in Stage 1 now
  inputTab = inputTab %>%
    mutate(Initial            = ifelse(!is.na(neutral) & neutral, 5, 3)) %>% ## ADDED "& neutral" to the line below - October 10,2024
    mutate(Rule_Initial       = ifelse(!is.na(neutral) & neutral, ifelse(datasets == "WHO", 1, 6), ifelse(datasets == "WHO", 5, 9))) %>%
    mutate(rule_initial_up    = (SOLO_SorR >= MIN_SOLO_COUNT_UPGRADE & !is.na(PPVc_SOLO_lb) & PPVc_SOLO_lb >= 0.25 & OR_SOLO > 1        & OR_SOLO_pval_FDR_sig)) %>%
    mutate(Initial      = ifelse(rule_initial_up, 1, Initial)) %>%
    mutate(Rule_Initial = ifelse(rule_initial_up, ifelse(datasets == "WHO", 2, 7), Rule_Initial)) %>%
    mutate(rule_pncA_down     = (gene == PZA_GENE  & !is.na(PPV_SOLO)     & PPV_SOLO < 0.4       & PPV_SOLO_ub < 0.75 & datasets == "WHO" & stage == 1)) %>%
    mutate(Rule_Initial = ifelse(rule_pncA_down, 3, Rule_Initial)) %>%
    mutate(rule_pncA_up       = (gene == PZA_GENE  & SOLO_R >= 2          & PPVc_SOLO >= PPVC_UPGRADE_THRESHOLD & stage == 1)) %>% ## Used to be PPV
    ## NOTE: mutation != LOF_LABEL was intentionally removed from the rule below
    mutate(rule_newGenes_up   = (paste(drug, gene, sep = "_") %in% paste(NEW_PAIRS$drug, NEW_PAIRS$gene, sep = "_") & SOLO_R >= 2 & PPVc_SOLO >= PPVC_UPGRADE_THRESHOLD & stage == 1)) %>%
    mutate(Initial      = ifelse(rule_pncA_down & Initial == 3, 4, ifelse((rule_pncA_up | rule_newGenes_up) & Initial == 3, 2, Initial))) %>%
    mutate(Rule_Initial = ifelse(rule_pncA_up | rule_newGenes_up, ifelse(datasets == "WHO", 4, 8), Rule_Initial)) 
  ## Reconcile initial confidence grading between the WHO and the ALL datasets by first separating them once again
  Tab0 = inputTab %>%
    dplyr::filter(datasets == "WHO")
  inputTab = inputTab %>%
    dplyr::filter(datasets == "ALL") %>%
    full_join(Tab0, by = c("variant", "drug", "gene", "mutation", "neutral", "stage"), suffix = c("_ALL", "_WHO")) %>%
    mutate(Initial_WHO = replace_na(Initial_WHO, 3)) ## Not all variants are initially WHO-graded!
  inputTab = inputTab %>%
    mutate(    Initial    = ifelse(Initial_WHO == Initial_ALL         , Initial_WHO , NA)) %>%
    mutate(    datasets   = ifelse(Initial_WHO == Initial_ALL         , ifelse(!is.na(neutral) & neutral, "WHO", "ALL+WHO"), NA)) %>%
    mutate(Initial  = ifelse(Initial_WHO == 3 & Initial_ALL != 3,                                 Initial_ALL , Initial ),
           datasets = ifelse(Initial_WHO == 3 & Initial_ALL != 3,                                 "ALL"       , datasets)) %>%
    mutate(Initial  = ifelse( is.na(Initial_ALL) | (Initial_ALL == 3 & Initial_WHO != 3),         Initial_WHO , Initial ),
           datasets = ifelse( is.na(Initial_ALL) | (Initial_ALL == 3 & Initial_WHO != 3),         "WHO"       , datasets)) %>%
    mutate(Initial  = ifelse( pmax(Initial_WHO, Initial_ALL) <= 2 & Initial_WHO != Initial_ALL,   Initial_WHO , Initial ),
           datasets = ifelse( pmax(Initial_WHO, Initial_ALL) <= 2 & Initial_WHO != Initial_ALL,   "WHO"       , datasets)) %>%
    mutate(Initial  = ifelse( Initial_WHO == 4 & Initial_ALL <= 2,                                INITIAL_FLAG, Initial ),
           datasets = ifelse( Initial_WHO == 4 & Initial_ALL <= 2,                                "FLAG"      , datasets))
  ## Extract additional grading criteria from PREV_VERSION, then prepare to compute the final grades and record additional grading criteria
  ## Note that only one rule is applied per variant-drug pair; those to which a rule has been applied are marked by setting anyRule to TRUE
  inputTab = inputTab %>%
    mutate(Final = Initial, `Additional grading criteria` = NA, Rule_Final = LAST_INITIAL_RULE, anyRule = FALSE)
  ## Upgrade to grade 4 variants initially graded 5 by set C only (i.e. literature only) or previous version guidance
  inputTab = inputTab %>%
    mutate(rule_literature       = ((setC_WHO | prev_version_WHO) & !(setA_WHO | setB_WHO | setD_WHO | setE_WHO) & Initial == 5)) %>%
    applyNamedRule("rule_literature")
  ## Downgrade to grade 2 variants initially graded 1 by the ALL dataset only
  inputTab = inputTab %>%
    mutate(rule_AllOnly     = (Initial_ALL == 1 & Initial_WHO == 3)) %>%
    applyNamedRule("rule_AllOnly")
  ## Mark those variants which had discrepant AwR grades in the two datasets as "Based on WHO dataset"
  ## Update (FEB 7, 2026): the indicator only applies when we have grade 2 dominating grade 1, not the other way around
  tempDescription = ifelse(is.na(inputTab$Initial_WHO) | is.na(inputTab$Initial_ALL) | inputTab$Initial_WHO != inputTab$Initial_ALL, WHO_BASED, NA)
  inputTab = inputTab %>%
    ## mutate(rule_WHOBased    = (pmax(Initial_WHO, Initial_ALL) <= 2 & Initial_WHO != Initial_ALL)) %>%
    mutate(rule_WHOBased    = (Initial_WHO == 2 & Initial_ALL == 1)) %>%
    applyNamedRule("rule_WHOBased", description = tempDescription)
  ## Mark those variants which triggered a flag at the initial classification as "Manual check required"
  inputTab = inputTab %>%
    mutate(rule_FLAG        = (Initial == INITIAL_FLAG)) %>%
    applyNamedRule("rule_FLAG")
  ## Upgrade rpoB borderline mutations to grade 1; note: this rule overrides any previously applied rules!
  tempDescription = ifelse(is.na(inputTab$Final) | inputTab$Final > 1, BORDERLINE, NA)
  inputTab = inputTab %>%
    mutate(rule_rpoB_borderline = (gene == RRDR_GENE & mutation %in% BORDERLINE_MUTATIONS)) %>%
    applyNamedRule("rule_rpoB_borderline", description = tempDescription)
  ## Downgrade any silent variant to grade 4, then set to NA their SOLO counts and related statistics!
  inputTab = inputTab %>%
    mutate(rule_silent = (effect_ALL %in% SILENT_EFFECTS                                                                    & Initial == 3)) %>%
    applyNamedRule("rule_silent") %>%
    mutate(across(matches("SOLO"), ~ifelse(Rule_Final == ruleEndNum("rule_silent"), NA, .)))
  ## Upgrade any non-silent variant in the RRDR region to grade 2
  inputTab = inputTab %>%
    mutate(rule_RRDR = (gene == RRDR_GENE & (pos1_ALL %in% RRDR_INT | pos2_ALL %in% RRDR_INT) & !effect_ALL %in% SILENT_EFFECTS & Initial == 3)) %>%
    applyNamedRule("rule_RRDR")
  ## If a pooled LoF is graded 1 or 2, then any grade 3 LoF mutation in the same drug-gene combination is upgraded to grade 2
  inputTab = inputTab %>%
    group_by(gene, drug) %>%
    mutate(rule_LoF_candidate = (any(mutation == LOF_LABEL & Final <= RESISTANCE_GRADE_MAX))) %>%
    ungroup() %>%
    mutate(rule_LoF = (rule_LoF_candidate & mutation != LOF_LABEL & effect_ALL %in% POOLED_EFFECTS[[LOF_LABEL]]                     & Initial == 3)) %>%
    applyNamedRule("rule_LoF")
  initRuleCrossRes = ruleEndNum("cross_res") - NUM_CROSS_RES_RULES + 1L
  ## Apply the cross-resistance rules the first time around
  inputTab = inputTab %>%
    applyCrossResistanceRules(iteration = 1, initRule = initRuleCrossRes)
  ## Upgrade to grade 2 any variant that is "recognized as DR marker" in a WHO-endorsed assay, as well as any variant with allelic exchange evidence
  assayTab = auxData$assayTab
  inputTab = inputTab %>%
    full_join(assayTab, by = c("drug", "variant")) %>%
    mutate(assay = convertToLogical(assay),   anyRule = convertToLogical(anyRule), rule_Assay = (assay & (is.na(Initial) | Initial == 3))) %>%
    applyNamedRule("rule_Assay")
  allelicTab = auxData$allelicTab
  inputTab = inputTab %>%
    full_join(allelicTab, by = c("drug", "variant")) %>%
    mutate(allelic = convertToLogical(allelic), anyRule = convertToLogical(anyRule), rule_Allelic = (allelic & (is.na(Initial) | Initial == 3)))
  tempDescription = ifelse(is.na(inputTab$Final) | inputTab$Final > 2, ALLELIC_EXP, NA)
  inputTab = inputTab %>%
    applyNamedRule("rule_Allelic", description = tempDescription) %>%
    select(-assay, -allelic)
  ## If any of the rules after the cross-resistance one has upgraded a pooled LOF variant to grade 2, we apply it one more time to the corresponding variants
  inputTab = inputTab %>%
    group_by(gene, drug) %>%
    mutate(rule_LoF_re_candidate = any(mutation == LOF_LABEL & Final <= RESISTANCE_GRADE_MAX & Rule_Final >= initRuleCrossRes)) %>%
    ungroup() %>%
    mutate(rule_LoF_re = (rule_LoF_re_candidate & mutation != LOF_LABEL & effect_ALL %in% POOLED_EFFECTS[[LOF_LABEL]]              & Initial == 3)) %>%
    applyNamedRule("rule_LoF_re")
  ## Apply the cross-resistance rules the second time around
  inputTab = inputTab %>%
    applyCrossResistanceRules(iteration = 2, initRule = initRuleCrossRes)
  ## Upgrade a specific variant to grade 2 based on Literature evidence (PMID 32571824)
  inputTab = inputTab %>%
    mutate(rule_PZA_lit = (variant == PZA_SPEC_VAR                                                                         & Initial == 3)) %>%
    applyNamedRule("rule_PZA_lit", description = describePMIDs(32571824))
  prevGuidanceRule = ruleEndNum("rule_prev_guidance")
  if (!any(inputTab$variant == DCS_SPEC_VAR)) {
    inputTab = inputTab %>%
      bind_rows(tibble(drug = "Cycloserine", variant = DCS_SPEC_VAR, anyRule = FALSE, Initial = 3, Final = 2, Rule_Final = prevGuidanceRule,
                       `Additional grading criteria` = describePMIDs(c(22912881, 27064254))))
  }
  ## Upgrade to grade 2 based on previous guidance to distinguish between "evidence in version 2" and "guidance before version 2" (e.g. Miotto ERJ2017)
  inputTab = inputTab %>%
    mutate(rule_prev_guidance = (!is.na(`Additional grading criteria`) & `Additional grading criteria` %in% PREV_GUIDANCE   & Initial == 3)) %>%
    applyNamedRule("rule_prev_guidance")
  extraTab = auxData$extraTab
  inputTab = inputTab %>%
    left_join(extraTab, by = c("drug", "variant")) %>%
    mutate(rule_prev_evidence = (!is.na(Final_prev_version) & Final_prev_version < 3 & (!effect_ALL %in% INFRAME_EFFECTS) & Initial == 3)) %>%
    applyExpertRule("rule_prev_evidence", description = PREV_EVIDENCE, finalGrade = FINAL_FLAG, finalRule = prevGuidanceRule + 0.5)
  ## Add comment column to a specified list of mutations or mutation categories, provided they were initially graded 3 and no other rule has applied:
  commentCategoryTab = auxData$commentCategoryTab
  commentLoF         = auxData$commentLoF
  commentSingleTab   = auxData$commentSingleTab
  inputTab = inputTab %>%
    full_join(commentLoF,         by = c("drug", "gene")) %>%
    mutate(comment = ifelse(!(effect_ALL %in% POOLED_EFFECTS[[LOF_LABEL]]), NA, comment)) %>%
    full_join(commentCategoryTab, by = c("drug", "gene", "Final"   ), suffix = c(".x", ".z")) %>%
    adjustDuplicateColumns(suffixes = c(".x", ".z"), add = TRUE) %>%
    full_join(commentSingleTab  , by = c("drug", "gene", "mutation"), suffix = c(".x", ".z")) %>%
    adjustDuplicateColumns(suffixes = c(".x", ".z"), add = TRUE) %>%
    dplyr::filter(!is.na(variant))
  ## Specify PMIDs that lead to a downgrade by the first 'proper' rule:
  inputTab = inputTab %>%
    mutate(`Additional grading criteria` = ifelse(drug == "Bedaquiline"      & Rule_Final == ruleEndNum("rule_literature"), describePMIDs(c(28031270, 34503982)), `Additional grading criteria`)) %>%
    mutate(`Additional grading criteria` = ifelse(drug %in% LEV_MXF          & Rule_Final == ruleEndNum("rule_literature"), describePMIDs(28137812),              `Additional grading criteria`)) %>%
    mutate(`Additional grading criteria` = ifelse(drug %in% FIRST_LINE_DRUGS & Rule_Final == ruleEndNum("rule_literature"), describePMIDs(32143680),              `Additional grading criteria`)) %>%
    mutate(`Additional grading criteria` = ifelse(drug == "Cycloserine"      & Rule_Final == ruleEndNum("rule_literature"), describePMIDs(c(27064254, 28971867)), `Additional grading criteria`))
  ## Remove variant addition variables
  inputTab = inputTab %>%
    select(-starts_with("add_"))
  inputTab = inputTab %>%
    mutate(temp = paste("from",   GRADES[Final_prev_version], "to", GRADES[Final])) %>%
    mutate(OverallFlag = ifelse(is.na(Final_prev_version), paste("New", GRADES[Final]),
                                ifelse(Final_prev_version == Final, "No change",
                                       ifelse((Final_prev_version <= RESISTANCE_GRADE_MAX & Final > 3) | (Final_prev_version > 3 & Final <= RESISTANCE_GRADE_MAX), paste("Switch", temp),
                                         ifelse(abs(Final_prev_version - 3) < abs(Final - 3), paste("Up", temp), paste("Down", temp)))))) %>%
    select(-temp)
  ## Lastly, convert the numerical grades to their description and save the file
  inputTab = inputTab %>%
    mutate(Initial = GRADES[Initial], Final = GRADES[Final])
  inputTab = inputTab %>%
    mutate(variant = ifelse(variant == "RRDR_LoF", "RRDR", variant))
  inputTab
}

#' Write graded variant results to a CSV file
#' @noRd
writeGradingOutput = function(inputTab, stage, LoF = TRUE, outDir = ".") {
  outFilename = paste(GRADED_CATALOGUE_PREFIX, Sys.Date(), sep = "_")
  outFilename = paste0(outFilename, ifelse(LoF, WITHLOFS_SUFFIX, ""), "_Stage", stage, ".csv")
  write_csv(inputTab, file = file.path(outDir, outFilename))
}

#' Format a list of PubMed IDs as an evidence string
#' @noRd
describePMIDs = function(PMIDs) {
  output = paste0("Literature evidence (PMID ", paste(PMIDs, collapse = "; "), ")")
  output = make_Evidence(output)
  output
}

#' Apply cross-resistance rules between drug pairs
#' @noRd
applyCrossResistanceRules = function(inputTab, iteration = 1, initRule = 0) {
  ## Set all numeric columns except for the identifying ones to NA
  numColumns = colnames(inputTab)[sapply(inputTab, is.numeric)]
  numColumns = setdiff(numColumns, c("stage", "Initial", "Final"))
  for (setting in 1:nrow(CROSS_RES)) {
    curSetting = CROSS_RES %>% 
      slice(setting)
    curDrugs = c(curSetting$drug1, curSetting$drug2)
    extra_variants = inputTab %>%
      group_by(variant) %>%
      mutate(add_variant = (any(sum(drug %in% curDrugs) == 1 & gene %in% curSetting$genes[[1]] & Final <= RESISTANCE_GRADE_MAX & !(effect_ALL %in% POOLED_EFFECTS[[LOF_LABEL]])))) %>%
      ungroup() %>%
      dplyr::filter(add_variant & drug %in% curDrugs) %>%
      mutate(drug = curDrugs[3 - match(drug, curDrugs)], Initial = 3, `Additional grading criteria` = NA_character_) %>%
      mutate(across(all_of(numColumns), ~NA_real_)) %>%
      mutate(across(starts_with('rule'), ~{ FALSE })) %>%
      mutate(anyRule = FALSE)
    inputTab = inputTab %>%
      bind_rows(extra_variants) %>%
      group_by(variant) %>%
      mutate(rule_cross_res = (gene %in% curSetting$genes[[1]] & any(drug %in% curDrugs & Final <= RESISTANCE_GRADE_MAX) & drug %in% curDrugs & Initial == 3)) %>%
      applyExpertRule("rule_cross_res", description = curSetting$description, finalGrade = 2, finalRule = initRule + (iteration - 1)/2) %>%
      ungroup()
    initRule = initRule + 1
  }
  inputTab
}

#' Apply a single expert grading rule to a variant table
#' @noRd
applyExpertRule = function(inputTab, ruleColumn, description = NA, finalGrade = NA, finalRule = NA, applyAlways = FALSE) {
  if (applyAlways) {
    inputTab = inputTab %>%
      mutate(applyRule = .data[[ruleColumn]])
  } else {
    inputTab = inputTab %>%
      mutate(applyRule = (!anyRule & .data[[ruleColumn]]))
  }
  if (!all(is.na(description))) { ### OCT 8, 2025: Change to allow description to be a vector; only use non-missing values when the rule applies
    inputTab = inputTab %>%
      mutate(`Additional grading criteria` = ifelse(applyRule & !is.na(description), description, `Additional grading criteria`))
  }
  if (!is.na(finalGrade)) {
    inputTab = inputTab %>%
      mutate(Final      = ifelse(applyRule, finalGrade, Final))
  }
  if (!is.na(finalRule)) {
    inputTab = inputTab %>%
      mutate(Rule_Final = ifelse(applyRule, finalRule, Rule_Final))
  }
  inputTab = inputTab %>%
    mutate(anyRule = or(anyRule, applyRule)) %>%
    select(-applyRule)
  inputTab
}
