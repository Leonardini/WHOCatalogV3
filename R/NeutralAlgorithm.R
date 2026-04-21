#' Identify neutral variants using a multi-set PPV-based algorithm
#' @param masterTab A data frame of genotype-phenotype data for all samples and drugs.
#' @param litTab A data frame of literature-based neutral mutations (loaded from literature_neutrals.csv).
#' @param catalogueFile Path to the previous-version WHO catalogue (.xlsx); used to identify resistance mutations for set B onwards.
#' @inheritParams mainDriver
#' @export
neutralAlgorithm = function(masterTab, litTab,
                            NON_DATABASE_DIRECTORY = system.file("extdata", package = "SOLOport"),
                            catalogueFile = file.path(NON_DATABASE_DIRECTORY, "who_catalogue_prev_version.xlsx")) {
  for (set in c("A", "B")) {
    if (set == "B") {
      ## Goal: remove the sample-drug pairs with a resistance mutation from the previous version, including based on the expert rules
      masterTab = masterTab %>%
        applyCatalogue(catalogueFile = catalogueFile, minMAF = MAF_THRESHOLD_REGULAR, lowMAFHet = TRUE,
                       minQ = QUALITY_THRESHOLD_STRICT, lowQHet = TRUE,
                       LoF = TRUE, version = PREV_VERSION)
      ## Computing the rows to be marked with rm = TRUE (Final <= 2); also converting rm to a logical variable (meaning NA -> FALSE)
      masterTab = masterTab %>%
        mutate(rm = (Final <= 2L)) %>%
        mutate_at("rm", convertToLogical)
      # Variants in some selected tier 1 genes that have an LoF mutation are also RMs
      masterTab %<>%
        mutate(rm = or(rm, paste0(drug, gene, sep = "_") %in% paste0(NEW_PAIRS_RM$drug, NEW_PAIRS_RM$gene, sep = "_") 
                       & effect %in% POOLED_EFFECTS[["LoF"]] & tier == 1))
      ## Mark any sample-drug pair with at least one of the variants being a resistance mutation
      masterTab %<>%
        group_by(sample_id, drug) %>%
        mutate(rm_in_sample = any(rm)) %>%
        ungroup()
    }
    ## Compute the PPVs
    masterTab %<>%
      computePPVs(removeRM = (set == "B"))
    ## Place the variant-drug combinations with PPV_ub < threshold into the set
    masterTab %<>%
      mutate(set = (!is.na(variant) & !(variant == "missing") & (!is.na(PPV_ub)) & (PPV_ub < PPV_UB_THRESHOLD)))
    ## Change the set name according to its letter
    colnames(masterTab)[ncol(masterTab)] = paste0("set", set)
    ## Extract the neutral variants into a separate file
    curNeutral = extractNeutral(masterTab, set)
    ## Remove generic auxiliary variables from the master table
    masterTab %<>%
      select(-PPV, -PPV_ub)
  }
  ## Goal: complete set C by adding literature mutations
  ## Merge with the master table
  masterTab %<>%
    full_join(litTab, by = c("drug", "variant"))
  ## Define and compute lit_mutation as a logical form of literature and remove literature
  masterTab %<>%
    mutate(lit_mutation = convertToLogical(literature)) %>%
    select(-literature)
  ## Define set C as the union of sets A, B and lit_mutations
  masterTab %<>%
    mutate(setC = (setA | setB | lit_mutation))
  ## Extract the neutral variants into a separate file
  curNeutral = extractNeutral(masterTab, "C")
  for (set in c("D", "E")) {
    ## filter out silent, tier 2 variants, and those in previous sets (C, D?), plus isolates containing a het
    subsetTab = masterTab %>%
      dplyr::filter(!(effect %in% SILENT_EFFECTS | tier == 2 | setC))
    if (set == "E") {
      subsetTab = subsetTab %>%
        dplyr::filter(!setD)
    }
    subsetTab %<>%
      removeSamplesWithHets()
    ## Compute the PPV_solo for this subset
    subsetTab %<>%
      computePPVs(removeRM = FALSE, solo = TRUE, restrict = TRUE)
    ## Add the computed values back to the master table 
    masterTab %<>%
      left_join(subsetTab, by = c("drug", "variant"))
    ## Place the variant-drug combinations with PPV_ub < threshold into the set
    masterTab %<>%
      mutate(set = ((!is.na(PPV_ub)) & (PPV_ub < PPV_UB_THRESHOLD)))
    ## Change the set name according to its letter
    colnames(masterTab)[ncol(masterTab)] = paste0("set", set)
    ## Extract the neutral variants into a separate file
    curNeutral = extractNeutral(masterTab, set)
    ## Remove generic auxiliary variables from the master table
    masterTab %<>%
      select(-PPV, -PPV_ub)
  }
  masterTab %<>%
    mutate(prev_version = FALSE)
  ## Convert the variable to a logical one
  masterTab %<>%
    mutate_at(c("prev_version", paste0("set", LETTERS[1:5]), "lit_mutation"), convertToLogical)
  ## Create the overall list of mutations
  masterTab %<>%
    mutate(setF = (setC | setD | setE | prev_version))
  ## Extract the neutral variants into a separate file
  curNeutral = extractNeutral(masterTab, "F")
  masterTab
}

#' Compute the positive predictive value (PPV) for each variant-drug pair
#' If removeRM = TRUE, the table is first filtered to samples without any resistance mutation.
#' If solo = TRUE, only variants that are solos within their sample-drug pair are kept.
#' If restrict = TRUE, returns only the filtered table with variants, drugs and PPV stats.
#' @noRd
computePPVs = function(inputTab, removeRM = TRUE, solo = FALSE, restrict = FALSE) {
  auxTab = inputTab
  if (removeRM && "rm_in_sample" %in% colnames(auxTab)) {
    auxTab %<>%
      dplyr::filter(!rm_in_sample)
  }
  if (solo) {
    ## Only keep the mutations that are solos within their respective sample-drug pair
    auxTab %<>%
      group_by(sample_id, drug) %>%
      mutate(N = n()) %>%
      ungroup() %>%
      mutate(solo = (N == 1)) %>%
      dplyr::filter(solo) %>%
      select(-solo, -N)
  }
  ## Compute the PPV and its confidence intervals for each variant-drug pair
  auxTab %<>%
    group_by(drug, variant) %>%
    mutate(total = sum(!het), totalR = sum((!het) & phenotype == "R")) %>%
    ## mutate(total = n(), totalR = sum(phenotype == "R")) %>%
    slice(1) %>%
    ungroup
  PPVTab = auxTab %>%
    distinct(totalR, total, .keep_all = FALSE) %>%
    mutate(testResult = map2(totalR, total, safeBinomTest))
  auxTab %<>%
    left_join(PPVTab, by = c("totalR", "total")) %>%
    mutate(PPV = map_dbl(testResult, ~{.$estimate}), PPV_ub = map_dbl(testResult, ~{.$conf.int[2]})) %>%
    select(drug, variant, PPV, PPV_ub)
  if (restrict) {
    return(auxTab)
  } else { ## Incorporate the computed values into the original input
    inputTab %<>%
      left_join(auxTab, by = c("drug", "variant"))
    return(inputTab)
  }
}

#' Extract neutral variant-drug combinations for a given neutral set
#' @noRd
extractNeutral = function(inputTab, setName) {
  ## Create the full variable name
  varName = paste0("set", setName)
  ## Extract the corresponding variable
  varMask = inputTab %>%
    select(one_of(varName)) %>%
    pull()
  ## Extract the variant-drug combinations in the set into a separate table
  neutralTab = inputTab %>%
    dplyr::filter(varMask) %>%
    group_by(drug, variant) %>%
    slice(1) %>%
    ungroup()
  ## Special case for set F: also keep track of all the other sets
  if (setName == "F") {
    neutralTab %<>%
      select(drug, variant, setA, setB, setC, setD, setE, prev_version, lit_mutation)
  } else {
    neutralTab %<>%
      select(drug, variant)
  }
  neutralTab
}

#' Write neutral mutation output files for all six sets (A through F)
#' NOTE: correctness relies on set marks (setA..setF, lit_mutation, prev_version) being frozen once
#' assigned in neutralAlgorithm — they must never be mutated afterwards.
#' @noRd
writeNeutralOutputs = function(masterTab, safe = FALSE, outDir = ".") {
  for (set in c("A", "B", "C", "D", "E", "F")) {
    write_csv(extractNeutral(masterTab, set), file.path(outDir, paste0("neutral_mutations_WHO_", set, ".csv")))
  }
  if (safe) {
    ## Extract the complete collection of marks for each variant-drug pair
    allSets = masterTab %>%
      group_by(drug, variant) %>%
      slice(1) %>%
      ungroup() %>%
      select(drug, variant, setA, setB, setC, lit_mutation, setD, setE, prev_version, setF)
    # Record the complete collection of marks
    write_csv(allSets, file.path(outDir, "allVariantsWithSetMarks.csv"))
  }
}
