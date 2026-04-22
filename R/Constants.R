## Index of the previous and current version
CURR_VERSION = 3L
PREV_VERSION = CURR_VERSION - 1L
EXTRACTION_ID = "2025-12-16T22-09-33.644305_jr_639ca1c8b567dc6a6e49b9097bd7c369c74bb0ba0a5749c1d9e886890abaa333"
PREFIX_STRING = "DATASET"
NUM_STAGES = 3L
LARGE_NUMBER = 2.5e6 ## Used for column type guessing

## List of all drugs and their short names 
DRUG_LIST   = c("Amikacin", "Bedaquiline", "Capreomycin", "Cycloserine", "Delamanid", "Ethambutol", 
                "Ethionamide", "Isoniazid", "Kanamycin", "Levofloxacin", "Linezolid", "Moxifloxacin", 
                "Para-Aminosalicylic Acid", "Pyrazinamide", "Rifampicin", "Streptomycin", "Clofazimine", "Pretomanid")
SHORT_NAMES = c("AMI", "BDQ", "CAP", "DCS", "DLM", "EMB", 
                "ETO", "INH", "KAN", "LEV", "LZD", "MXF", 
                "PAS", "PZA", "RIF", "STM", "CFZ", "PMD")
NUM_DRUGS = length(DRUG_LIST)
EXTRA_DRUGS = c("Ofloxacin", "Prothionamide")
EXTRA_NAMES = c("OFX", "PTO")

## Special-status drugs
LEV_MXF          = DRUG_LIST[match(c("LEV", "MXF")                     , SHORT_NAMES)]
INH_ETO          = DRUG_LIST[match(c("INH", "ETO")                     , SHORT_NAMES)]
BDQ_CFZ          = DRUG_LIST[match(c("BDQ", "CFZ")                     , SHORT_NAMES)]
DLM_PMD          = DRUG_LIST[match(c("DLM", "PMD")                     , SHORT_NAMES)]
NUM_CROSS_RES_RULES = 4 ## This needs to match the number of groups computed just above!
FIRST_LINE_DRUGS = DRUG_LIST[match(c("EMB", "INH", "PZA", "RIF")       , SHORT_NAMES)]
GROUP_3_DRUGS    = DRUG_LIST[match(c("BDQ", "CFZ", "DLM", "INH", "PZA"), SHORT_NAMES)]

## Special-status genes
BDQ_GENE     = "Rv0678"
CFZ_GENE     = "pepQ"
BDQ_CFZ_GENE = c(BDQ_GENE, CFZ_GENE)
CAP_GENE     = "tlyA"
DLM_GENE     = c("ddn", "fbiA", "fbiB", "fbiC", "fgd1", "Rv2983")
DLM_PMD_GENE = DLM_GENE
ETO_GENE     = c("ethA", "mshA", "mshC")
LEV_MXF_GENE = c("gyrA", "gyrB")
INH_GENE     = "katG"
INH_ETO_GENE = c("inhA","fabG1")
PZA_GENE     = "pncA"
RRDR_GENE    = "rpoB"
STM_GENE     = "gid"
PAS_GENE     = "thyA"
DCS_GENE     = "cycA"

COMPENSATORY = paste0("ahpC_c.-", c("48G>A", "51G>A", "52C>A", "52C>T", "54C>T", "57C>T", "72C>T", "76T>A", "77T>A", "81C>T"))

NEW_PAIRS    = tibble(drug = DRUG_LIST[match(c("INH", "BDQ", "CFZ", rep(c("DLM", "PMD"), length(DLM_GENE)), "PZA"), SHORT_NAMES)], 
                      gene = c(INH_GENE, rep(c(BDQ_GENE, DLM_GENE), each = 2), PZA_GENE))
NEW_PAIRS_RM = rbind(NEW_PAIRS, 
                     tibble(drug = DRUG_LIST[match("PAS", SHORT_NAMES)], gene = PAS_GENE), 
                     tibble(drug = DRUG_LIST[match("CAP", SHORT_NAMES)], gene = CAP_GENE),
                     tibble(drug = DRUG_LIST[match("STM", SHORT_NAMES)], gene = STM_GENE),
                     tibble(drug = DRUG_LIST[match("DCS", SHORT_NAMES)], gene = DCS_GENE),
                     tibble(drug = DRUG_LIST[match("ETO", SHORT_NAMES)], gene = ETO_GENE))

## Special-status gene sets
ADDITIONAL_GENES   = list(DLM = DLM_GENE, PMD = DLM_GENE, INH = INH_GENE, PZA = PZA_GENE)
EXTENDED_ADD_GENES = c(ADDITIONAL_GENES, list(BDQ = BDQ_CFZ_GENE, CFZ = BDQ_CFZ_GENE, CAP = CAP_GENE, ETO = ETO_GENE[1], STM = STM_GENE))
ADDITIONAL_GENES   = c(ADDITIONAL_GENES, list(BDQ = BDQ_GENE, CFZ = BDQ_GENE, PAS = PAS_GENE))
SUB_EXTENDED_GENES = EXTENDED_ADD_GENES[c("ETO", "INH", "PZA", "STM")]
EXCLUDE_BDQ_GENES  = c("atpE",  "pepQ")
SPECIAL_BDQ_GENES  = c("mmpL5")
STRATIFY_BDQ_GENES = c(SPECIAL_BDQ_GENES, "mmpS5")

## Special-status variants
BORDERLINE_MUTATIONS = c("p.Leu430Pro", "p.Asp435Tyr", "p.His445Leu", "p.His445Asn", "p.His445Ser", "p.Leu452Pro", "p.Ile491Phe")
## INFLATED_PPV_VARS    = c("p.Asn98fs"  , "p.Cys46Arg" , "p.Cys46fs"  , "p.Gln51fs"  , "p.Ile67Ser" , "p.Leu142fs" , "p.Met146Thr", "p.Pro48fs")
PZA_SPEC_VAR         = paste(PZA_GENE, "p.Ile31Thr" , sep = "_")
DCS_SPEC_VAR         = paste(DCS_GENE, "p.Gly122Ser", sep = "_")
SPECIAL_DLM_VAR      = paste(DLM_GENE[1], "p.Leu49Pro", sep = "_")

## Special-status effects
LOF_LABEL            = "LoF"
POOLED_EFFECTS       = list(LoF = c(paste0(c("feature", "transcript"), "_ablation"), "frameshift", "start_lost", "stop_gained"))
INFRAME_EFFECTS      = paste0("inframe_", c("deletion", "insertion"))
ADDITIONAL_EFFECTS   = c(POOLED_EFFECTS[["LoF"]], INFRAME_EFFECTS, "missense_variant")
BDQ_EFFECTS          = c(ADDITIONAL_EFFECTS, "stop_lost")
SILENT_EFFECTS          = c("initiator_codon_variant", "stop_retained_variant", "synonymous_variant")
NEUTRAL_LOGICAL_COLUMNS = c(paste0("set", LETTERS[1:5]), "lit_mutation", "prev_version")
UPSTREAM_VAR         = "upstream_gene_variant"
MISSING_VARIANT      = "missing"
INDEL_EFFECTS        = c(INFRAME_EFFECTS, paste0(c("feature", "transcript"), "_ablation"), "frameshift")

## Useful constants and thresholds
RRDR_INT              = 426:452
PPV_UB_THRESHOLD      = 0.1
SIG_THRESHOLD         = 0.05
FDR_THRESHOLD         = 0.05
CI_COEFFICIENT        = 1.96

MAF_THRESHOLD_STRICT      = 0.90
MAF_THRESHOLD_REGULAR     = 0.75
MAF_THRESHOLD_RELAXED     = 0.25
QUALITY_THRESHOLD_STRICT  = 1000
QUALITY_THRESHOLD_REGULAR = 500
QUALITY_THRESHOLD_RELAXED = 250
MAF_THRESHOLDS     = c(strict = MAF_THRESHOLD_STRICT,     regular = MAF_THRESHOLD_REGULAR,     relaxed = MAF_THRESHOLD_RELAXED)
QUALITY_THRESHOLDS = c(strict = QUALITY_THRESHOLD_STRICT, regular = QUALITY_THRESHOLD_REGULAR, relaxed = QUALITY_THRESHOLD_RELAXED)

RESISTANCE_GRADE_MAX   = 2L
MIN_SOLO_COUNT_UPGRADE = 5L
PPVC_UPGRADE_THRESHOLD = 0.5

## Grades used to prioritise variants
GRADES = paste0(c("", "not "), "assoc w ") %>% 
  str_to_sentence() %>%
  outer(c("R", "R - Interim"), function(x, y) { paste0(x, y) }) %>%
  c("Uncertain significance", "Manual check") %>%
  magrittr::extract(c(1, 3, 5, 4, 2, 6)) %>%
  enframe(value = "grading") %>%
  mutate(grading = paste0(name, ") ", grading)) %>%
  pull(grading)
INITIAL_FLAG     = 6
FINAL_FLAG       = 6
MAX_GRADE        = 6
LAST_INITIAL_RULE = 15L
FIRST_RULE_NUM    = LAST_INITIAL_RULE + 1L

## Lineage groupings
LINEAGE_OTHER = "Other"
LINEAGE_ANY   = "Any"

## Output file names
WITHLOFS_SUFFIX               = "_withLoFs"
STATS_FILE_PREFIX             = "Stats_"
GRADED_CATALOGUE_PREFIX       = "Final_graded_algorithm_catalogue"
COMPLETE_DATASET_FILENAME     = "CompleteDataset.csv"
COMPLETE_DATASET_WHO_FILENAME = "CompleteDatasetWHO.csv"

## Process-specific tables
## PHENO_GROUPS    = tibble(category_phenotype = c("ALL", "WHO", "CC", "CC-ATU"), group = c("MAIN", "MAIN", "CC", "ATU"))
PHENO_GROUPS       = tibble(category_phenotype = c("ALL", "WHO"), group = "MAIN")
BAD_VAR_DRUG_PAIRS = tibble(drug = c("Isoniazid", "Rifampicin"), variant = c("katG_p.Ser315Thr", "rpoB_p.Ser450Leu"))
CONVERTED_FILES    = paste0(c(paste0("neutral_mutations_catalogue_v", PREV_VERSION), "new_variant_matched_to_old", 
                                     paste0("v", PREV_VERSION, "_grades")), ".csv") ## USED TO BE v1 in both filenames

#' Prefix a string with "Rule:"
#' @noRd
make_Rule = function(x) { x = paste("Rule:", x); x }
#' Prefix a string with "Evidence:"
#' @noRd
make_Evidence = function(x) { x = paste("Evidence:", x); x }
#' Expand drug short names to full names in a string
#' @noRd
spellOut = function(x) { for (ind in 1:NUM_DRUGS) { x = str_replace(x, SHORT_NAMES[[ind]], DRUG_LIST[[ind]]) } ; x }

## Additional information supporting the old grades
OLD_LOF    = "Indel or premature stop codon (LoF)"
ALL_ONLY   = make_Evidence("Evidence from ALL dataset only")
PASS_2     = "Algorithm pass 2"
RRDR       = make_Rule("RRDR")
BORDERLINE = make_Evidence("Borderline")
PROBLEMATIC_CRITERION = make_Evidence("WHO-endorsed gDST assay")
IGNORED_CRITERIA = c(OLD_LOF, ALL_ONLY, PASS_2, RRDR, BORDERLINE)

## Additional information supporting the new grades
LIT_EVIDENCE = make_Evidence("Literature")
## INFLATION     = "Potentially inflated PPV"
## SET_C_ONLY    = "Neutrals defined by setC (literature) only"
## PREV_NEUTRALS = "previous WHO guidance" ## USED TO BE V1_NEUTRALS!
WHO_BASED     = make_Evidence("WHO grading takes priority")
MANUAL_CHECK  = "Manual check required"
SILENT        = make_Rule("Silent mutation")
LoF_MUT       = make_Rule("LoF")
WHO_ASSAY     = make_Evidence("Recognized as DR marker through WHO-endorsed assay")
ALLELIC_EXP   = make_Evidence("Selection")
PREV_GUIDANCE = "Previous WHO guidance"
PREV_EVIDENCE = make_Evidence(paste("Catalogue version", PREV_VERSION)) ## USED TO BE V1_EVIDENCE!

CROSS_RES = bind_rows(tibble(description = "Fluoroquinolones", drug1 = LEV_MXF[1], drug2 = LEV_MXF[2], genes = list(LEV_MXF_GENE)),
                      tibble(description = "INH-ETO",          drug1 = INH_ETO[1], drug2 = INH_ETO[2], genes = list(INH_ETO_GENE)),
                      tibble(description = "BDQ-CFZ",          drug1 = BDQ_CFZ[1], drug2 = BDQ_CFZ[2], genes = list(BDQ_CFZ_GENE)),
                      tibble(description = "DLM-PMD",          drug1 = DLM_PMD[1], drug2 = DLM_PMD[2], genes = list(DLM_PMD_GENE))) %>%
  mutate(description = make_Evidence(spellOut(paste(description, "cross-resistance"))))

## Constants used for SOLO algorithm results
MAX_ITER = Inf
UCODE    = 0L
RCODE    = 1L
SCODE    = -1L
CODE_KEY = c("S", "U", "R") %>% magrittr::set_names(c(SCODE, UCODE, RCODE))
HET_TAB  = tibble(class = c("S", "S", "U", "U", "R"), het = c(FALSE, TRUE, FALSE, TRUE, FALSE))
HET_CNT  = table(HET_TAB$class)
## and the final tabulation
STAGES   = c(1, 2, 3, 12, 123)

## Epistasis-related variants
EXCLUDE_SET = vector("list", 5)
EXCLUDE_SET[["AMI"]]     =                         paste0("eis_c.", "-14C>T")
EXCLUDE_SET[["AMI_EXT"]] =   EXCLUDE_SET[["AMI"]]
EXCLUDE_SET[["KAN"]]     = c(EXCLUDE_SET[["AMI"]], paste0("eis_c.", "-10G>A"))
EXCLUDE_SET[["KAN_EXT"]] = c(EXCLUDE_SET[["KAN"]], paste0("eis_c.", c("-12C>T", "-37G>T", "-8delC")))
EXCLUDE_SET[["BOTH"]]    =                         paste0("rrs_n.", c("1401A>G", "1402C>T", "1484G>T"))

