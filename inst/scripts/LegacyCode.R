library(exact2x2)

## Auxiliary function for converting all Stata input files into CSV format, keeping names but changing extensions
convertAllStataFiles = function(inputDir = "SOLO Algorithm Input files/STATA DTA files/") {
  LF = list.files(inputDir, pattern = ".dta$", full.names = TRUE)
  for (fpath in LF) {
    Tab = haven::read_dta(fpath)
    print(basename(fpath))
    write_csv(Tab, stringr::str_replace(fpath, ".dta$", ".csv"))
  }
}

oldCrossCheck = function() {
  Tab1 = read_csv("genome_indices_from Tim_22Apr2021.csv", guess_max = 2^20)
  Tab2 = read_csv("UNITAID_FIND_WHO_SimplifiedTable_FINAL_31Mar2021.csv", guess_max = 2^20)
  Tab3 = read_csv("Genome_Indices_FINAL_22Apr2021.csv", guess_max = 2^20)
  coln = colnames(Tab2)
  firstRow = Tab2 %>% 
    slice(1)
  coln[str_sub(coln, 1, 1) == "X"] = firstRow[str_sub(coln, 1, 1) == "X"]
  colnames(Tab2) = unlist(coln)
  Tab2 = Tab2 %>%
    dplyr::filter(c(FALSE, rep(TRUE, nrow(Tab2) - 1)))
  total2 = Tab2 %>%
    select(starts_with("MUT")) %>%
    mutate_all(as.integer) %>%
    rowSums(na.rm = TRUE)
  Tab2 = Tab2 %>%
    dplyr::filter(total2 > 0)
  Tab2 = Tab2 %>%
    mutate(variant = `variant (common_name)` %>% str_remove_all("\\ .*"))
  Tab2E = Tab2 %>%
    inner_join(Tab1, by = c("variant" = "variant", "drug" = "drug"))
  Tab2E = Tab2E %>%
    select(c(drug:`Genome position`, `FINAL CONFIDENCE GRADING`:found_in_S)) %>%
    select(-c(tier, `variant (common_name)`)) %>%
    select(gene, variant, `Genome position`, ref, alt, found_in_R, found_in_S, `FINAL CONFIDENCE GRADING`, everything())
  Tab2E = Tab2E %>%
    rename(gene_name = gene, genome_index = `Genome position`, ref_nt = ref, alt_nt = alt, Total_R = found_in_R, 
           Total_S = found_in_S, Conf_Grade = `FINAL CONFIDENCE GRADING`)
  Tab2E = Tab2E %>%
    mutate(mutation = str_sub(variant, nchar(gene_name) + 2)) %>%
    mutate(codon_number = str_extract(mutation, "^[A-Z]{1}[0-9]+[A-Z\\!]{1}$") %>% str_sub(2, -2) %>% as.integer()) %>%
    mutate(ref_aa = str_extract(mutation, "^[A-Z]"), alt_aa = str_extract(mutation, "[A-Z\\!]$")) %>%
    select(-mutation)
  Tab2E = Tab2E %>%
    select(gene_name, variant, codon_number, genome_index, ref_nt, alt_nt, ref_aa, alt_aa, Total_R, Total_S, everything())
  Tab2E = Tab2E %>%
    dplyr::filter(Conf_Grade != "combo")
  Tab2F = Tab2E %>%
    group_by(variant) %>%
    group_split() %>%
    lapply(function(x) {pivot_wider(x, names_from = drug, values_from = c(Total_R, Total_S, Conf_Grade))})
  Tab2G = do.call(bind_rows, Tab2F)
  TotalR = Tab2G %>%
    select(starts_with("Total_R")) %>%
    rowSums(na.rm = TRUE)
  TotalS = Tab2G %>%
    select(starts_with("Total_S")) %>%
    rowSums(na.rm = TRUE)
  Tab2G = Tab2G %>%
    mutate(Total_R = TotalR, Total_S = TotalS) %>%
    select(gene_name:alt_aa, Total_R, Total_S, Total_R_AMI:Conf_Grade_RIF, genome_index1:alt3)
  Tab2G = Tab2G %>%
    mutate(alt_codon_number = variant %>% str_extract("\\_[0-9]+\\_(ins|del)\\_") %>% str_remove_all("(ins|del)") %>% 
             str_remove_all("\\_") %>% as.integer() %>% add(2) %>% divide_by_int(3) %>% str_c("(", ., ")"))
  Tab3E = Tab3 %>%
    dplyr::filter(is.na(ref_aa) | is.na(alt_aa) | ref_aa != alt_aa)
  Vars3 = unique(Tab3E$variant)
  Vars2 = unique(Tab2G$variant)
  missingVars = setdiff(Vars2, Vars3)
  print(paste("There are", length(missingVars), "variants that only appear in the initial catalogue"))
  print(paste("They are:",  paste(missingVars, collapse = ", ")))
  omittedVars = setdiff(Vars3, Vars2)
  print(paste("There are", length(omittedVars), "variants that only appear in the final catalogue"))
  print(paste("They are:",  paste(omittedVars, collapse = ", ")))
  Tab2H = Tab2G %>%
    dplyr::filter(variant %in% Vars3) %>%
    arrange(variant, ref_nt, alt_nt)
  ### Source: https://mycobrowser.epfl.ch/releases (release 4, tab-separated files)
  extraTab1 = read_tsv("Mycobacterium_tuberculosis_H37Rv_txt_v4.txt", guess_max = 2^20) %>%
    select(c(Locus, Name)) %>%
    set_colnames(c("gene_locus", "gene_name"))
  ### Source: https://www.ddbj.nig.ac.jp/ddbj/code-e.html
  extraTab2 = read_tsv("AminoAcids.txt") %>%
    set_colnames(c("Abbreviation", "abbr", "name")) %>%
    select(c(Abbreviation, abbr)) %>%
    slice(-27) %>%
    bind_rows(tibble(Abbreviation = "Ter", abbr = "!"))
  Tab2H = Tab2H %>%
    left_join(extraTab1, by = "gene_name") %>%
    mutate_at(c("ref_aa", "alt_aa"), ~{MM = match(., extraTab2$abbr); out = extraTab2[MM, "Abbreviation"]; out})
  Tab3E = Tab3E %>%
    arrange(variant, ref_nt, alt_nt)
  coln = colnames(Tab2H)
  coln[11:55] = paste0(str_sub(coln[11:55], -3, -1), "_", str_sub(coln[11:55], 1,-5))
  coln[11:55] = str_remove(coln[11:55], "_Total")
  coln[56:64] = paste0("detail_", coln[56:64])
  colnames(Tab2H) = coln
  for (ind in 1:nrow(Tab2H)) {
    if (ind %% 1000 == 0) { print(ind) }
    curRow = Tab2H %>%
      slice(ind)
    curInds = curRow %>%
      select(c(detail_genome_index1, detail_genome_index2, detail_genome_index3))
    curInds = curInds[!is.na(curInds)]
    if (length(curInds) > 0) {
      Tab2H$genome_index[ind] = paste(curInds, collapse = "") %>% 
        as.numeric()
    }
    curCodons = curRow %>%
      select(c(codon_number, alt_codon_number))
    curCodons = curCodons[!is.na(curCodons)]
    if (length(curCodons) > 0) {
      stopifnot(length(curCodons) == 1)
      Tab2H$codon_number[ind] = paste(curCodons, collapse = "")
    }
  }
  Tab2H = Tab2H %>%
    select(-alt_codon_number)
  coln = colnames(Tab2H)
  stopifnot(colnames(Tab3E)[61] == "detail_ref2")
  Tab3F = Tab3E %>%
    dplyr::filter(!(variant %in% Tab2H$variant)) %>%
    group_by(variant) %>%
    mutate(N = n()) %>%
    ungroup %>%
    arrange(-N)
  Tab3H = Tab3E %>%
    dplyr::filter(variant %in% Tab2H$variant)
  ### Special exception for Rv1258c_c-23t
  Tab3H = Tab3H %>%
    dplyr::filter(variant != "Rv1258c_c-23t" | Total_S != 0)
  for (Col in coln) {
    print(Col)
    curColTrain = Tab3H %>%
      select(one_of(Col))
    curColTest = Tab2H %>%
      select(one_of(Col))
    if (Col == "genome_index") {
      stopifnot(all.equal.numeric(curColTrain, curColTest))
    } else {
      n1 = sum(curColTrain != curColTest, na.rm = TRUE)
    }
    if (n1 > 0) { stop(paste("There are", n1, "discrepancies")) }
    n2 = sum(is.na(curColTrain) & !is.na(curColTest))
    if (n2 > 0) { stop(paste("There are", n2, "missing values in training set only")) }
    n3 = sum(!is.na(curColTrain) & is.na(curColTest))
    if (n3 > 0) { print(paste("There are", n3, "missing values in testing set only")) }
  }
}

processTabs = function(inputFiles = c("Full_dataset_ALL.csv", "Full_dataset_WHO.csv")) {
  outputFiles = vector("list", length(inputFiles))
  for (ind in 1:length(inputFiles)) {
    inputFile = inputFiles[ind]
    print(paste("Processing input file", ind, ":", inputFile))
    Tab = read_csv(inputFile) %>%
      select(c("drug", "tier", "variant", starts_with("...")))
    coln = colnames(Tab)
    colnT = c(coln[1:3], Tab %>% 
                slice(1) %>% 
                select(-(1:3))) 
    colnT %<>% 
      unlist %>% 
      magrittr::set_names(NULL) %>%
      str_remove_all("MUT") %>% 
      str_remove_all("present as") %>% 
      str_remove_all("pheno") %>% 
      str_remove_all(" ")
    TabT = Tab %>%
      set_colnames(colnT) %>%
      slice(-1) %>%
      rowid_to_column("origIndex")
    mutInfo = TabT %>%
      select(c(origIndex, drug, variant))
    TabT %<>%
      select(-c(drug, variant)) %>%
      mutate_all(as.integer) %>%
      inner_join(mutInfo, by = "origIndex") %>%
      mutate(SOLO_S = SOLO_SorR - SOLO_R)
    TabT %<>%
      computeCatalogueStats()
    outputFile = str_replace(inputFile, ".csv", paste0("_Processed.csv"))
    write_csv(TabT, outputFile)
    outputFiles[[ind]] = TabT
  }
  outputFiles
}

## Compute the set of mutations that are in set A but not in set B
postprocessNeutralLists = function() {
  masterTab = read_csv("allVariantsWithSetMarks.csv")
  microTab = masterTab %>% 
    dplyr::filter(setA & !setB) %>%
    group_by(variant, drug) %>% 
    slice(1) %>%
    select(variant, drug)
  ## Write that set of mutations into a separate file
  write_csv(microTab, "neutral_mutations_WHO_A_not_B.csv")
}

safeBinomMeld = function(x, y, z, w) {
  if (y == 0 || w == 0) {
    output = list(estimate = NA, conf.int = c(NA, NA))
  } else {
    output = binomMeld.test(x1 = x, n1 = y, x2 = z, n2 = w, parmtype = "ratio", alternative = "greater")
  }
  output
}

oldCatalogueStats = function(TabT) {
  TabT %<>%
    mutate(      LRP = pmap(list(present_S, present_S + absent_S, present_R, present_R + absent_R), safeBinomMeld)) %>%
    mutate( LRP_SOLO = pmap(list(   SOLO_S,    SOLO_S + absent_S,    SOLO_R,    SOLO_R + absent_R), safeBinomMeld)) %>%
    mutate(      LRN = pmap(list( absent_S, present_S + absent_S,  absent_R, present_R + absent_R), safeBinomMeld)) %>%
    mutate( LRN_SOLO = pmap(list( absent_S,    SOLO_S + absent_S,  absent_R,    SOLO_R + absent_R), safeBinomMeld))
}

postprocessGradedStats = function(fast = FALSE, correct_all = TRUE) {
  matchTab = read_csv("new_variant_matched_to_old.csv", guess_max = LARGE_NUMBER, show_col_types = FALSE)
  for (name in unique(PHENO_GROUPS$group)) {
    if (!fast || name == "MAIN") {
      if (name == "MAIN") {
        Tab0 = read_csv("Graded_Stats_WHO.csv" , guess_max = LARGE_NUMBER, show_col_types = FALSE)
        Tab1 = read_csv("Graded_Stats_MAIN.csv", guess_max = LARGE_NUMBER, show_col_types = FALSE)
        Tab01 = full_join(Tab0, Tab1, by = c("drug", "variant")) %>%
          mutate(useWHO = or(useWHO.x, useWHO.y), neutral = or(neutral.x, neutral.y)) %>%
          select(-useWHO.x, -useWHO.y, -neutral.x, -neutral.y)
        Tab01 %<>%
          mutate_at("useWHO", ~{na_if(., FALSE)}) %>%
          group_by(drug, variant) %>%
          mutate(useWHO = any(useWHO)) %>%
          ungroup() %>%
          mutate(useWHO = useWHO | neutral | (!is.na(Initial.x) & Initial.x != 3 & Initial.y == 3)) %>%
          mutate_at("useWHO",
                    ~{ifelse(!is.na(Initial.y) & (Initial.x == Initial.y | Initial.y != 3) & is.na(.), FALSE, .)}) %>%
          mutate(special_case = is.na(useWHO)) %>%
          mutate_at("useWHO", ~{replace_na(., FALSE)}) %>%
          select(drug, variant, useWHO, special_case) %>%
          group_by(drug, variant) %>%
          slice(1) %>%
          ungroup()
        useWHO = Tab01 %>% 
          dplyr::filter(useWHO)
        useALL = Tab01 %>%
          dplyr::filter(!useWHO)
        finalTab = bind_rows(Tab0 %>% select(-useWHO) %>% inner_join(useWHO, by = c("drug", "variant")), 
                             Tab1 %>% select(-useWHO) %>% inner_join(useALL, by = c("drug", "variant")))
        finalTab %<>%
          mutate(datasets = ifelse(useWHO, "WHO", "ALL")) %>%
          mutate(downgrade = (special_case & Final == 1 & datasets == "ALL")) %>%
          mutate_at("Additional grading criteria_v1", ~{ifelse(downgrade, "Downgraded to interim based on WHO dataset", .)}) %>%
          mutate_at("Final"                         , ~{ifelse(downgrade, 2,                                            .)}) %>%
          select(-special_case, -downgrade) %>%
          mutate_at("Additional grading criteria_v1", ~{ifelse(Final %in% c(1, 5) & datasets == "ALL" & is.na(.), ALL_ONLY, .)}) %>%
          mutate_at("Final",                       ~{ifelse(`Additional grading criteria_v1` == ALL_ONLY, 2, .)})
      } else {
        finalTab = read_csv(paste0("Graded_Stats_", name, ".csv"), guess_max = LARGE_NUMBER, show_col_types = FALSE)
      }
      finalTab %<>%
        mutate_at("Additional grading criteria_v1", ~{paste(., ifelse(algorithm_pass == 2 & Final %in% c(1,5), "Algorithm pass 2", ""))}) %>%
        mutate_at("Final",                       ~{ifelse(Final == 1 & algorithm_pass == 2, 2, .)}) %>%
        mutate(Initial_Confidence_Grading = GRADES[Initial], Final_Confidence_Grading = GRADES[Final]) %>%
        select(-Initial, -Final)
      finalTab %<>%
        rename("variant_category_v2" = variant) %>%
        left_join(matchTab, by = c("drug", "variant_category_v2")) %>%
        mutate(Present_in_Catalogue_v1 = !is.na(variant_category_v1)) %>%
        rename("variant" = variant_category_v2, "variant_v1_nomenclature" = variant_category_v1)
      write_csv(finalTab, paste0(paste("Final_graded_algorithm_catalogue", name, Sys.Date(), "Leonid", 
                                       "Correct", ifelse(correct_all, "All", "SOLO"), sep = "_"), ".csv"))
    }
  }
}

EXTRA_LOF = FALSE # Set to TRUE if you want to generate information for the other two types of pooled indels as well!
if (EXTRA_LOF) {
  POOLED_EFFECTS[["inframe"]] = paste0("inframe_", c("insertion", "deletion"))
  POOLED_EFFECTS[["LoF_ALL"]] = c(POOLED_EFFECTS[["LoF"]], POOLED_EFFECTS[["inframe"]])
}

GROUP_123_DRUGS  = DRUG_LIST[match(c("BDQ", "INH", "PZA")              , SHORT_NAMES)]

if (stage == 2 & relax == FALSE & version == 2 & geno == "ALL") {
  miniSubset = curSubset %>%
    filter(selected & phenotype == "R" & drug == "CFZ") %>%
    distinct(sample_id, .keep_all = TRUE)
  write_csv(miniSubset, "Group2_75_RIFALL_CFZ_TP_samples.csv")
}

geneMap = fullDataset %>% select(gene, drug) %>% distinct()

## CODE TO CHECK FOR THE EXISTENCE OF PAIRS OF VARIANTS WITH ONE HET AND ONE NON-HET AFFECTING THE SAME POSITION IN THE SAME GENE
fullDataset %<>%
  bind_rows(orphanData) %>%
  mutate(het_relaxed = (`max(af)` < MAF_THRESHOLD_RELAXED | variant == "missing"), het_strict = (`max(af)` < MAF_THRESHOLD_STRICT | variant == "missing")) %>%
  select(sample_id, drug, variant, phenotype, gene, mutation, effect, pos1, pos2, `max(af)`, het, het_relaxed, het_strict)
miniDataset = fullDataset %>%
  group_by(sample_id, drug, gene, pos1) %>%
  mutate(N = n()) %>%
  filter(N > 1) %>%
  ungroup()
miniDataset %<>%
  group_by(sample_id, drug, pos1) %>%
  mutate(K = n_distinct(het)) %>%
  filter(K > 1 & !is.na(pos1)) %>%
  ungroup() %>%
  arrange(sample_id, drug, pos1)
altDataset = fullDataset %>%
  mutate(nuc_change = str_detect(mutation, ">")) %>%
  group_by(sample_id, drug, gene, pos1, nuc_change) %>%
  mutate(S = sum(`max(af)`, na.rm = TRUE)) %>%
  filter(S > 1) %>%
  ungroup() %>%
  arrange(sample_id, drug, pos1)

## Sanity check: most of the isolates containing a URM should have a resistant phenotype
## miniTab = table(masterTab %>% slice(1) %>% pull(phenotype), masterTab %>% slice(1) %>% pull(urm))
## print("Displaying the table of URM presence/absence by phenotype; one unit is a sample-drug combination")
## print(miniTab)

## From the NewCrossCheck function:
# resultsJoint %<>%
#   testAndRemove("Sens", "Sensitivity")
# resultsJoint %<>%
#   testAndRemove("sens_all_lb", "Sens_lb")
# resultsJoint %<>%
#   testAndRemove("sens_all_ub", "Sens_ub")
# resultsJoint %<>%
#   testAndRemove("Spec", "Specificity")
# resultsJoint %<>%
#   testAndRemove("spec_all_lb", "Spec_lb")
# resultsJoint %<>%
#   testAndRemove("spec_all_ub", "Spec_ub")

### Obsolete function from NeutralAlgorithm.R to convert version 1 to version 2 variants
convertVersions = function(safe = FALSE, NON_DATABASE_DIRECTORY = "STATA DTA files", DATA_DIRECTORY = "DATABASE EXTRACTION files", 
                           EXTRACTION_ID = NULL) {
  ## Check if the files already exist
  if (!safe & all(file.exists(CONVERTED_FILES))) {
    print("All the required files have already been created; exiting!")
    return()
  }
  ## Save the starting directory to move back there at the end and go into the matching directory
  initDir = getwd()
  setwd(paste(DATA_DIRECTORY, EXTRACTION_ID, "prev_version_matching/", sep="/"))
  ## Preprocess the first file by renaming its columns and creating a new one for compatibility with current version
  matchTab = read_csv(list.files()[1], guess_max = LARGE_NUMBER, show_col_types = FALSE) %>%
    rename(variant_category_prev_version = 'description', gene = 'resolved_symbol', mutation = 'variant_category', effect = 'predicted_effect') %>%
    mutate(variant_category_curr_version = paste(gene, mutation, sep = "_")) %>%
    select(-gene, -mutation)
  ## Keep only the non-silent mutations
  matchNonSilentTab = matchTab %>%
    dplyr::filter(!(effect %in% SILENT_EFFECTS))
  ## Extract revised list of neutral mutations
  oldTab = read_csv_2headers(paste0(NON_DATABASE_DIRECTORY, "/List_of_neutral_mutations_from_cat_ver_", PREV_VERSION, "_rev20Jan2023.csv"), prefix = FALSE)
  ## Rename the drugs according to their full names
  oldTab %<>%
    mutate_at("drug", ~{DRUG_LIST[match(., SHORT_NAMES)]})
  ## Remove any entries that were marked as not usable
  oldTab %<>%
    dplyr::filter(!str_starts(`FOR WHO CATALOGUE UPDATE VER. 2`, "NOT")) %>%
    select(drug, variant) %>%
    rename(variant_category_prev_version = "variant")
  ## Find the common substrate between those and the matched non-silent variants
  oldTab %<>%
    inner_join(matchNonSilentTab, by = "variant_category_prev_version") %>%
    rename(variant = "variant_category_curr_version")
  ## Read in the final report from the publication (NB: the csv file contains only the Mutation_catalogue sheet)
  altTab = read_csv_2headers(paste0(NON_DATABASE_DIRECTORY, "/WHO-UCN-GTB-PCI-2021.7-eng.csv"), prefix = FALSE) %>%
    mutate(across(starts_with("Present") | starts_with("Absent"), as.integer))
  ## Remove combo graded entries
  altTab %<>%
    dplyr::filter(`FINAL CONFIDENCE GRADING` != "combo") %>%
    rename(variant = `variant (common_name)`)
  ## Rename the drugs according to their full names
  altTab %<>%
    mutate_at("drug", ~{DRUG_LIST[match(., SHORT_NAMES)]})
  ## Remove the Genome position column, normalize variant names, remove empty ones, and rename some columns
  altTab %<>%
    select(-`Genome position`) %>%
    mutate_at("variant", ~{str_split_fixed(., " ", n = 2)[,1]}) %>%
    dplyr::filter(variant != "") %>%
    rename(variant_category_prev_version = "variant", final_grading_prev_version = `FINAL CONFIDENCE GRADING`) %>%
    rename(`Additional grading criteria_prev_version` = `Additional grading criteria`)
  ## Replace any problematic grading criteria with an empty string
  altTab %<>%
    mutate_at("Additional grading criteria_prev_version", ~{na_if(., PROBLEMATIC_CRITERION)})
  ## Check that all of the matched variants are present in the original table
  stopifnot(all(matchNonSilentTab$variant_category_prev_version %in% altTab$variant_category_prev_version))
  ## Merge the two tables; variant_category_prev_version is NA for the original variants that do not have a matched variant
  altTab %<>%
    full_join(matchNonSilentTab, by = "variant_category_prev_version", relationship = "many-to-many")
  ## Among the entries that do not have a matched variant, exclude all except those that we want to postprocess
  ## First, drop variants in the fprA gene as it is no longer on the list
  altTab %<>%
    dplyr::filter(!(is.na(variant_category_curr_version) & str_starts(variant_category_prev_version, "fprA")))
  ## Mark the variants that are neither insertions nor deletions
  altTab %<>%
    mutate(non_indel = !str_detect(variant_category_prev_version, "del") & !str_detect(variant_category_prev_version, "ins"))
  ## Drop the variants in certain genes that are neither insertions nor deletions, and have no matched variant
  for (gene in c("rpsL", "embR", "rrs", "embC", "ubiA")) {
    altTab %<>%
      dplyr::filter(!(is.na(variant_category_curr_version) & str_starts(variant_category_prev_version, gene) & non_indel))
  }
  ## Drop all other promoter and indel variants
  altTab %<>%
    dplyr::filter(!(is.na(variant_category_curr_version) & (!non_indel | str_detect(variant_category_prev_version, "-"))))
  ## Remove the non_indel column
  altTab %<>%
    select(-non_indel)
  ## Drop all the uncertain variants
  altTab %<>%
    dplyr::filter(!(is.na(variant_category_curr_version) & str_detect(variant_category_prev_version, "[MV]1[A-Z]")))
  ## Drop the gid_V110G and inhA_T4I variants
  altTab %<>%
    dplyr::filter(!(is.na(variant_category_curr_version) & variant_category_prev_version %in% c("gid_V110G", "inhA_T4I")))
  ## Manually recode the remaining pncA mutations
  xTab = tibble(variant_category_prev_version = paste0("pncA_"  , c("H71D",     "I31T",     "L116R",     "T135N")), 
                variant_category_curr_version = paste0("pncA_p.", c("His71Asp", "Ile31Thr", "Leu116Arg", "Thr135Asn")))
  altTab %<>%
    recodeValues(xTab)
  ## Double-check that every entry is annotated with a variant in the current version now ### TEMPORARY FIX ADDED: MAY 27, 2025 
  stopifnot((!is.na(altTab$variant_category_curr_version)) | (str_starts(altTab$variant_category_prev_version, 'ndh_')))
  ## Keep a single representative per drug-current version variant combination
  ## New addition: remove any variants in current version that correspond to more than one variant in previous version
  altTab %<>%
    group_by(drug, variant_category_curr_version) %>%
    mutate(N = n()) %>%
    filter(N == 1) %>%
    ungroup() %>%
    select(-N)
  ## Recode the final grading values for certain specific variants
  yTab = tibble(variant_category_curr_version = c("ethA_p.Met1?", "katG_p.Met1?", "pncA_p.Met1?"), 
                final_grading_prev_version = GRADES[c(1, 3, 1)])
  altTab %<>%
    recodeValues(yTab)
  ## Save a version with old and new variants and previous version's confidence grading
  matchingTab = altTab %>% 
    select(drug, variant_category_prev_version, variant_category_curr_version, final_grading_prev_version, `Additional grading criteria_prev_version`)
  ## Then extract additional data from the corresponding additional dataset
  extraTab = read_csv(paste0(NON_DATABASE_DIRECTORY, "/assay_mutations_24Apr2023.csv"), guess_max = LARGE_NUMBER, show_col_types = FALSE) %>%
    mutate_at("drug", ~{DRUG_LIST[match(., SHORT_NAMES)]}) %>%
    rename("variant_category_curr_version" = "variant") %>%
    mutate(`Additional grading criteria_curr_version` = PROBLEMATIC_CRITERION)
  matchingTab %<>%
    full_join(extraTab, by = c("variant_category_curr_version", "drug")) %>%
    mutate_at("Additional grading criteria_prev_version", ~{ifelse(`Additional grading criteria_curr_version` == PROBLEMATIC_CRITERION, 
                                                                   PROBLEMATIC_CRITERION, .)})
  ## Save another version with new variants and additional grading criteria
  gradingTab = altTab %>% 
    rename(variant = "variant_category_curr_version") %>% 
    select(drug, `Additional grading criteria_prev_version`, final_grading_prev_version, variant)
  ## Go back to the starting directory and write the files
  setwd(initDir)
  ## Save the results into appropriate files
  write_csv(oldTab,      file = CONVERTED_FILES[1])
  write_csv(matchingTab, file = CONVERTED_FILES[2])
  write_csv(gradingTab,  file = CONVERTED_FILES[3])
}

## Auxiliary function for processing files (usually derived from Excel spreadsheets) with two header rows
## If prefix = TRUE, prepends the first header row (propagated rightwards) to the second header row
## Otherwise, directly uses the second header row whenever it is non-missing
read_csv_2headers = function(inputFile, prefix = FALSE) {
  initTab = read_csv(inputFile, guess_max = LARGE_NUMBER, show_col_types = FALSE, name_repair = "unique_quiet")
  ## Merge the second row of headers into the first one
  initHeaders = initTab %>% 
    slice(1)
  coln = colnames(initTab)
  if (prefix) { 
    for (ind in 2:length(coln)) {
      if (is.na(coln[ind])) {
        coln[ind] = coln[ind - 1]
      }
    }
    coln[!is.na(initHeaders)] = paste(coln[!is.na(initHeaders)], initHeaders[!is.na(initHeaders)], sep = "_")
  } else {
    coln[!is.na(initHeaders)] = initHeaders[!is.na(initHeaders)]
  }
  initTab %<>%
    slice(-1) %>%
    set_colnames(coln)
  initTab
}

## Auxiliary function for testing that all the groups defined by groupingVars agree on consistentVars in a Table
testConsistentOld = function(Table, groupingVars, consistentVars) {
  L = length(consistentVars)
  if (L == 0) { return(TRUE) }
  cnames = colnames(Table)
  stopifnot(all(groupingVars %in% cnames))
  stopifnot(all(consistentVars %in% cnames))
  groupedTab = Table %>%
    group_by(across(all_of(groupingVars)))
  checks   = rep(FALSE, L) %>%
    magrittr::set_names(consistentVars)
  problems = vector("list", L) %>%
    magrittr::set_names(consistentVars)
  for (ind in 1:L) {
    curVar = consistentVars[ind]
    miniTab = groupedTab %>%
      select(any_of(curVar)) %>%
      distinct() %>%
      mutate(N = n()) %>%
      ungroup()
    if (!(all(miniTab$N == 1))) {
      problems[[ind]] = miniTab %>%
        dplyr::filter(N > 1)
    } else {
      checks[ind] = TRUE
    }
  }
  output = list(checks = checks, problems = problems)
}

## Obsolete pieces from SensSpec.R
# } else {
#   for (curName in names(fullTab)) {
#     curTab = fullTab[[curName]] %>%
#       mutate(mut_R = has_MOI & phenotype == "R", mut_S = has_MOI & phenotype == "S") %>%
#       slice(1) %>%
#       ungroup() %>%
#       mutate_at(c("mut_R", "mut_S"), sum) %>%
#       slice(1) %>%
#       select(mut_R:mut_S) %>%
#       mutate(drug = curDrug, group = curName, .before = 1) %>%
#       mutate(relevant_MOI = NA)
#     miniTab %<>%
#       bind_rows(curTab)
#   }

## COMMENTED OUT 13 FEB 2025
# if (version == CURR_VERSION && relaxed) {
#   gradedTab %<>% 
#     rename(Final_Confidence_Grading = `FINAL CONFIDENCE GRADING`)
# }

## Obsolete code for v1; commented out on 21 July 2025
# } else {
#   fullDataset %<>%
#     mutate(LOF_candidate = ((drug_short == "CAP" & gene == CAP_GENE & effect %in% setdiff(POOLED_EFFECTS[["LoF"]], "start_lost")) |
#                               (gene %in% unlist(SUB_EXTENDED_GENES) & effect %in% c(POOLED_EFFECTS[["LoF"]], INFRAME_EFFECTS))))
# }

## Obsolete piece from ComputeStatsNew.R - original version of the FDR correction, with a more explicit alternative
## mutate(OR_p_rank      = rank(OR_p)     , OR_p_max      = max(OR_p_rank     [OR_p     /OR_p_rank      <= FDR_threshold/k])) %>%
## mutate(OR_SOLO_p_rank = rank(OR_SOLO_p), OR_SOLO_p_max = max(OR_SOLO_p_rank[OR_SOLO_p/OR_SOLO_p_rank <= FDR_threshold/k])) %>%
## mutate(     OR_pvalue_threshold =      OR_pvalue[     OR_pvalue_rank ==      OR_pvalue_max][1]) %>%
## mutate(OR_SOLO_pvalue_threshold = OR_SOLO_pvalue[OR_SOLO_pvalue_rank == OR_SOLO_pvalue_max][1]) %>%

## Obsolete pieces from NeutralAlgorithm.R
## Newly added: variants in NEW_PAIRS that have a coding, non-silent mutation are also URMs - commented out on July 21, 2025
# SILENT_OR_NONCODING_EFFECTS = c(SILENT_EFFECTS, UPSTREAM_VAR, "non_coding_transcript_exon_variant", "missing")
# masterTab %<>%
#   mutate(urm = or(urm, paste(drug, gene, sep = "_") %in% paste(NEW_PAIRS_URM$drug, NEW_PAIRS_URM$gene, sep = "_") & 
#                     (!(effect %in% SILENT_OR_NONCODING_EFFECTS))))
## The code below is now obsolete - commented out on July 21, 2025
# ## Variants in BDQ_GENE that have an effect among the BDQ_EFFECTS are also URMs
# masterTab %<>%
#   mutate(urm = or(urm, gene == BDQ_GENE & effect %in% BDQ_EFFECTS))

fname = paste(NON_DATABASE_DIRECTORY, "urmLatest.csv", sep="/")
if (!file.exists(fname)) {
  ## If the file does not exist, try to create it from the previous extraction
  print("The URM file does not exist; trying to create it from the previous extraction")
  result = prepPreviouslyRAndS(NON_DATABASE_DIRECTORY = NON_DATABASE_DIRECTORY)
}
urmTab = read_csv(fname, show_col_types = FALSE, guess_max = LARGE_NUMBER)
## NEWLY ADDED: Variants that are LoF mutations are URMs too, but must be treated in a special way due to disaggregation in v3+
LoFTab = urmTab %>%
  filter(str_detect(variant, "LoF")) ## COMPLETE THIS! NOTE THAT THEY NEED TO BE APPLIED TO PMD AS WELL!
## MAKE ALL THE DRUG-GENE PAIRS USED EXPLICIT THROUGHOUT THE CODE!
## Apply the grading rule to these: if a pooled LoF is graded 1 or 2, any LoF mutation in the same drug-gene combination is a URM
urmTab %<>%
  filter(!str_detect(variant, "LoF")) %>%
  mutate(urm = TRUE)
## Merge with the master table
masterTab %<>%
  left_join(urmTab, by = c("drug", "variant"))

# The code below is now obsolete
# ## Import previous version mutations, convert them to current version, then mark as neutral
# convertVersions(safe = safe, NON_DATABASE_DIRECTORY = NON_DATABASE_DIRECTORY, DATA_DIRECTORY = DATA_DIRECTORY, 
#                 EXTRACTION_ID = EXTRACTION_ID)

## Obsolete definitions from Constants.R
## EXTRA_GENES  = c("mshA", "panD", "Rv1979c") ## OBSOLETE!
## EXTRA_GENES   = c("ald", "mshA", "mshC")
## SPECIAL_RM_GENES = c(CAP_GENE, CFZ_GENE, ETO_GENE[1], PAS_GENE, STM_GENE, EXTRA_GENES, unlist(ADDITIONAL_GENES)) %>% magrittr::set_names(NULL)

## Alternative based on a bedaquiline issue: if skipBDQfromSA is TRUE, remove the Bedaquiline entries corresponding to a list of sample_id's
if (skipBDQfromSA) {
  curDir = getwd()
  setwd(DATA_DIRECTORY)
  badIDs = read_csv("query_result_2023-04-13T15_28_45.986784Z.csv", show_col_types = FALSE, guess_max = LARGE_NUMBER) %>%
    pull(`Sample ID`)
  for (name in names(fullDataset)) {
    fullDataset[[name]] %<>%
      filter(!(sample_id %in% badIDs & drug == "Bedaquiline"))
  }
  setwd(curDir)
}

## Block at the end of the mainDriver function that adds orphan mutations to the final output file - now obsolete
finalResult = prepareOutputForPaolo(gradedFilename, LoF = LoF)
setwd(paste(str_remove(DATA_DIRECTORY, "/$"), str_remove(EXTRACTION_ID, "/$"), "orphan_genotypes/", sep = "/"))
orphanTab = read_csv(list.files()[1], guess_max = LARGE_NUMBER, show_col_types = FALSE) %>%
  rename(gene = 'resolved_symbol', mutation = 'variant_category', effect = 'predicted_effect', drug = 'drug_name') %>%
  mutate(variant = paste(gene, mutation, sep = "_"))
if (orphanByDrug) {
  orphanTab %<>%
    select(variant, drug, sample_id) %>%
    group_by(variant, drug)
} else {
  orphanTab %<>%
    select(variant, sample_id) %>%
    group_by(variant)
}
orphanTab %<>%
  mutate(Present_NoPheno = n_distinct(sample_id)) %>%
  slice(1) %>%
  ungroup() %>%
  select(-sample_id)
print(paste("Adding in the", nrow(orphanTab), "orphan mutations"))
setwd(OUTPUT_DIRECTORY)
if (orphanByDrug) {
  finalResult %<>%
    full_join(orphanTab, by = c("variant", "drug"))
} else {
  finalResult %<>%
    full_join(orphanTab, by = "variant")
}
finalResult %<>% mutate_at("Present_NoPheno", ~{replace_na(., 0)})
write_csv(finalResult, paste0(gradedFilename, "_withOrphans", ifelse(orphanByDrug, "_Stratified", ""), ".csv"))

## Obsolete option from prepMask: aS = TRUE masks algorithmic S variants (classified S by the SOLO algorithm)
if (aS) {
  inputTab %<>%
    dplyr::filter(is.na(class) | class != "S")
}

prepareOutputForPaolo = function(inputFile, LoF = FALSE) {
  Tabs = vector("list", 2)
  Names = c("WHO", "ALL")
  inputTab = read_csv(inputFile, show_col_types = FALSE, guess_max = LARGE_NUMBER) %>%
    rename(Datasets = datasets, Neutral_masked = neutral, OR_SOLO_FE_sig = OR_SOLO_pval_FDR_sig) %>%
    rename(PPV_conditional_SOLO = PPVc_SOLO, PPV_conditional_SOLO_lb = PPVc_SOLO_lb, PPV_conditional_SOLO_ub = PPVc_SOLO_ub) %>%
    rename(Absent_R = absent_R, Absent_S = absent_S, Present_R = present_R, Present_S = present_S, Present_SOLO_R = SOLO_R, Present_SOLO_SR = SOLO_SorR) %>%
    rename(Additionalgradingcriteria = `Additional grading criteria_prev_version`, Final_Confidence_Grading = Final, Initial_Confidence_Grading = Initial) %>%
    select(-correctAll, -correctSOLO, -FDR_threshold, -OR_pvalue, -OR_pval_FDR_sig, -OR_pval_max, -OR_pval_rank, -OR_SOLO_pval_max, -stage) %>%
    select(-starts_with("set"))
  for (index in 1:2) {
    curTab = inputTab %>%
      filter(Datasets == Names[index])
    curColn = colnames(curTab)
    colnames(curTab) = paste0(Names[index], "_", curColn)
    Tabs[[index]] = curTab
  }
  finalTab = full_join(Tabs[[1]], Tabs[[2]], by = c("WHO_drug" = "ALL_drug", "WHO_variant" = "ALL_variant"))
  write_csv(finalTab, paste0("Merged_Graded_Stats_", Sys.Date(), ifelse(LoF, "_withLoFs", ""), ".csv"))
  finalTab
}

## This function prepares the list of mutations previously associated with resistance as well as those not previously associated with it
prepPreviouslyRAndS = function(inputFile = "WHO-UCN-TB-2023.7-eng.xlsx", NON_DATABASE_DIRECTORY = "STATA DTA files") {
  initTab = read_excel_2headers(paste(NON_DATABASE_DIRECTORY, inputFile, sep = "/"), postfix = TRUE)
  coln = unique(colnames(initTab))
  initTab = initTab[, coln] %>%
    as_tibble()
  prevListR = initTab %>%
    filter(`FINAL CONFIDENCE GRADING` %in% GRADES[1:2]) %>%
    select(drug, variant)
  prevListS = initTab %>%
    filter(`FINAL CONFIDENCE GRADING` == GRADES[5] & !(effect %in% INDEL_EFFECTS)) %>%
    select(drug, variant)
  write_csv(prevListR, file = paste(NON_DATABASE_DIRECTORY, "rmLatest.csv",    sep = "/"))
  write_csv(prevListS, file = paste(NON_DATABASE_DIRECTORY, CONVERTED_FILES[1], sep = "/"))
  finalTab = initTab %>%
    mutate(final_grading_prev_version = `FINAL CONFIDENCE GRADING`, `Additional grading criteria_prev_version` = `Additional grading criteria applied`) %>%
    select(drug, `Additional grading criteria_prev_version`, final_grading_prev_version, variant)
  write_csv(finalTab,  file = paste(NON_DATABASE_DIRECTORY, CONVERTED_FILES[3], sep = "/"))
  output = list(R = prevListR, S = prevListS)
  output
}

applyCrossResistanceRules = function(inputTab, iteration = 1, initRule = 0) {
  ## Set all numeric columns except for the identifying ones to NA
  numColumns = colnames(inputTab)[sapply(inputTab, is.numeric)]
  numColumns = setdiff(numColumns, c("stage", "Initial", "Final"))
  ## Upgrade to grade 2 based on LEV-MFX cross-resistance (LEV-MFX on gyrA/gyrB, bilateral); first, add grade 1/2 non-LoF variants appearing in one drug
  extra_FQ_variants = inputTab %>%
    group_by(variant) %>%
    mutate(add_FQ_variant = (any(sum(drug %in% LEV_MXF) == 1 & gene %in% LEV_MXF_GENE & Final < 3 & !(effect_ALL %in% POOLED_EFFECTS[["LoF"]])))) %>%
    ungroup() %>%
    filter(add_FQ_variant & drug %in% LEV_MXF) %>%
    mutate(drug = LEV_MXF[3 - match(drug, LEV_MXF)], Initial = 3, `Additional grading criteria` = NA_character_) %>%
    mutate_at(numColumns, ~{ NA_real_ }) %>%
    mutate(across(starts_with('rule'), ~{ FALSE })) %>%
    mutate(anyRule = FALSE)
  inputTab = inputTab %>%
    bind_rows(extra_FQ_variants) %>%
    group_by(variant) %>%
    mutate(rule_FQ_cross_res = (gene %in% LEV_MXF_GENE  & any(drug %in% LEV_MXF & Final < 3) & drug %in% LEV_MXF           & Initial == 3)) %>%
    applyExpertRule("rule_FQ_cross_res",     description = FQ_CROSS_RES, finalGrade = 2, finalRule = initRule + (iteration - 1)/2) %>%
    ungroup()
  initRule = initRule + 1
  ## Upgrade to grade 2 based on INH-ETO cross-resistance (INH-ETO on inhA/fabG1, bilateral); first, add grade 1/2 non-LoF variants appearing in one drug
  extra_IE_variants = inputTab %>%
    group_by(variant) %>%
    mutate(add_IE_variant = (any(sum(drug %in% INH_ETO) == 1 & gene %in% INH_ETO_GENE & Final < 3 & !(effect_ALL %in% POOLED_EFFECTS[["LoF"]])))) %>%
    ungroup() %>%
    filter(add_IE_variant & drug %in% INH_ETO) %>%
    mutate(drug = INH_ETO[3 - match(drug, INH_ETO)], Initial = 3, `Additional grading criteria` = NA_character_) %>%
    mutate_at(numColumns, ~{ NA_real_ }) %>%
    mutate(across(starts_with('rule'), ~{ FALSE })) %>%
    mutate(anyRule = FALSE)
  inputTab = inputTab %>%
    bind_rows(extra_IE_variants) %>%
    group_by(variant) %>%
    mutate(rule_IE_cross_res = (gene %in% INH_ETO_GENE & any(drug %in% INH_ETO & Final < 3) & drug %in% INH_ETO          & Initial == 3)) %>%
    applyExpertRule("rule_IE_cross_res",      description = IE_CROSS_RES, finalGrade = 2, finalRule = initRule + (iteration - 1)/2) %>%
    ungroup()
  initRule = initRule + 1
  ## Upgrade to grade 2 based on BDQ-CFZ cross-resistance (BDQ-CFZ on Rv0678/pepQ, bilateral); first, add grade 1/2 non-LoF variants appearing in one drug
  ## CHANGED FROM UNILATERAL TO BILATERAL, 5 DECEMBER 2024
  extra_BC_variants = inputTab %>%
    group_by(variant) %>%
    mutate(add_BC_variant = (any(sum(drug %in% BDQ_CFZ) == 1 & gene %in% BDQ_CFZ_GENE & Final < 3 & !(effect_ALL %in% POOLED_EFFECTS[["LoF"]])))) %>%
    ungroup() %>%
    filter(add_BC_variant & drug %in% BDQ_CFZ) %>%
    mutate(drug = BDQ_CFZ[3 - match(drug, BDQ_CFZ)], Initial = 3, `Additional grading criteria` = NA_character_) %>%
    mutate_at(numColumns, ~{ NA_real_ }) %>%
    mutate(across(starts_with('rule'), ~{ FALSE })) %>%
    mutate(anyRule = FALSE)
  inputTab = inputTab %>%
    bind_rows(extra_BC_variants) %>%
    group_by(variant) %>%
    mutate(rule_BC_cross_res = (gene %in% BDQ_CFZ_GENE & any(drug %in% BDQ_CFZ & Final < 3) & drug %in% BDQ_CFZ        & Initial == 3)) %>%
    applyExpertRule("rule_BC_cross_res",      description = BC_CROSS_RES, finalGrade = 2, finalRule = initRule + (iteration - 1)/2) %>%
    ungroup()
  initRule = initRule + 1
  ## Upgrade to grade 2 based on DLM-PMD cross-resistance (DLM-PMD on ddn, fbiA, fbiB, fbiC, fgd1, Rv2983, bilateral); first, add grade 1/2 non-LoF variants appearing in one drug
  ## ADDED ON 5 DECEMBER 2024 [LC]
  extra_DP_variants = inputTab %>%
    group_by(variant) %>%
    mutate(add_DP_variant = (any(sum(drug %in% DLM_PMD) == 1 & gene %in% DLM_PMD_GENE & Final < 3 & !(effect_ALL %in% POOLED_EFFECTS[["LoF"]])))) %>%
    ungroup() %>%
    filter(add_DP_variant & drug %in% DLM_PMD) %>%
    mutate(drug = DLM_PMD[3 - match(drug, DLM_PMD)], Initial = 3, `Additional grading criteria` = NA_character_) %>%
    mutate_at(numColumns, ~{ NA_real_ }) %>%
    mutate(across(starts_with('rule'), ~{ FALSE })) %>%
    mutate(anyRule = FALSE)
  inputTab = inputTab %>%
    bind_rows(extra_DP_variants) %>%
    group_by(variant) %>%
    mutate(rule_DP_cross_res = (gene %in% DLM_PMD_GENE & any(drug %in% DLM_PMD & Final < 3) & drug %in% DLM_PMD       & Initial == 3)) %>%
    applyExpertRule("rule_DP_cross_res",      description = DP_CROSS_RES, finalGrade = 2, finalRule = initRule + (iteration - 1)/2) %>%
    ungroup()
  inputTab
}


## ---------------------------------------------------------------------------
## Legacy code from SensSpec.R / computeSensSpec (removed during D-section cleanup)

## Alternative epistasis exclusion: mark epistasis candidates as het before
## computing group min (requires excludeEpi_Regular/Relaxed computed in
## computeEpistasisStats; disabled pending full v3 epistasis implementation)
# if (version == CURR_VERSION) {
#   fullDataset = fullDataset %>%
#     mutate_at("het",         ~{ifelse(excludeEpi_Regular, TRUE, .)}) %>%
#     mutate_at("het_relaxed", ~{ifelse(excludeEpi_Relaxed, TRUE, .)})
# }

## Alternative extendedCandidate logic based on ADDITIONAL_EFFECTS / ADDITIONAL_GENES
## (replaced by the katG-specific INH version active in computeSensSpec)
# fullDataset = fullDataset %>%
#   mutate(extendedCandidate = (effect %in% ADDITIONAL_EFFECTS & (gene %in% unlist(ADDITIONAL_GENES[-1]) | (drug_short == "DLM" & gene %in% EXTENDED_ADD_GENES[["DLM"]])))) %>%
#   mutate_at("Group_Regular", ~{ifelse(RIF_geno_Regular & any(!het         & Final == 3 & extendedCandidate & . == 4), 3, .)}) %>%
#   mutate_at("Group_Relaxed", ~{ifelse(RIF_geno_Relaxed & any(!het_relaxed & Final == 3 & extendedCandidate & . == 4), 3, .)}) %>%
#   select(-extendedCandidate)
