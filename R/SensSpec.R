#' Load and prepare the complete dataset for sensitivity/specificity analysis
#' @noRd
loadSensSpecData = function(DATA_DIRECTORY, version, useOrphanData = TRUE, useLineageData = TRUE, useSublineageData = FALSE) {
  fullDataset = read_csv(paste0("Results/", EXTRACTION_ID, "/CompleteDataset.csv"), guess_max = LARGE_NUMBER, show_col_types = FALSE)
  if (useOrphanData) {
    orphanData = getOrphanData(DATA_DIRECTORY, EXTRACTION_ID)
    fullDataset %<>%
      bind_rows(orphanData)
  }
  if (useLineageData) {
    lineageData = getLineageData(DATA_DIRECTORY, EXTRACTION_ID, useSublineageData)
    fullDataset = fullDataset %>%
      left_join(lineageData)
  }
  fullDataset
}

#' Assign an epistasis group label to each sample based on a specified gene
#' @noRd
assignEpiGroup = function(df, geneName) {
  df %>%
    mutate(Group = ifelse( any(gene == geneName &  effect %in% POOLED_EFFECTS[["LoF"]] & !het_strict),                      "A",
                           ifelse(!any(gene == geneName & !(effect %in% c(SILENT_EFFECTS, UPSTREAM_VAR)) & Final %in% 1:3), "B", NA)))
}

#' Summarize epistasis group assignments across samples
#' @noRd
summarizeEpiGroup = function(df, groupByCols, extraCols) {
  df %>%
    dplyr::filter(!is.na(Group)) %>%
    group_by(across(all_of(c(groupByCols, extraCols)))) %>%
    summarise(N = n(), .groups = "drop")
}

#' Compute epistasis groupings and PPV tables for AMI, KAN, and BDQ
#' @noRd
computeEpistasisStats = function(fullDataset) {
  fullDataset %<>%
    mutate(excludeEpi_Candidate = FALSE, excludeEpi_Regular = FALSE, excludeEpi_Relaxed = FALSE) %>%
    mutate_at("excludeEpi_Candidate", ~{ifelse(drug_short == "AMI"            , any(!het_strict & gene == "eis"   & effect %in% POOLED_EFFECTS[["LoF"]]), .)}) %>%
    mutate_at("excludeEpi_Candidate", ~{ifelse(drug_short == "KAN"            , any(!het_strict & gene == "eis"   & effect %in% POOLED_EFFECTS[["LoF"]]), .)}) %>%
    mutate_at("excludeEpi_Candidate", ~{ifelse(drug_short %in% c("BDQ", "CFZ"), any(!het_strict & gene == "mmpL5" & effect %in% POOLED_EFFECTS[["LoF"]]), .)}) %>%
    mutate_at("excludeEpi_Regular"  , ~{ifelse(excludeEpi_Candidate & drug_short == "AMI"             & !het         & variant %in% EXCLUDE_SET[["AMI"]],  TRUE,
                                        ifelse(excludeEpi_Candidate & drug_short == "KAN"             & !het         & variant %in% EXCLUDE_SET[["KAN"]],  TRUE,
                                        ifelse(excludeEpi_Candidate & drug_short %in% c("BDQ", "CFZ") & !het         & gene    %in% BDQ_GENE & Final <= 2, TRUE, .)))}) %>%
    mutate_at("excludeEpi_Relaxed"  , ~{ifelse(excludeEpi_Candidate & drug_short == "AMI"             & !het_relaxed & variant %in% EXCLUDE_SET[["AMI"]],  TRUE,
                                        ifelse(excludeEpi_Candidate & drug_short == "KAN"             & !het_relaxed & variant %in% EXCLUDE_SET[["KAN"]],  TRUE,
                                        ifelse(excludeEpi_Candidate & drug_short %in% c("BDQ", "CFZ") & !het_relaxed & gene    %in% BDQ_GENE & Final <= 2, TRUE, .)))})
  fullDataset %<>% ungroup()
  epiTabs = vector("list", 4) %>%
    magrittr::set_names(c("AMI", "KAN", paste0("BDQ_", STRATIFY_BDQ_GENES)))
  for (drugName in c("AMI", "KAN")) {
    curEpiTab = fullDataset %>%
      dplyr::filter(drug_short == drugName) %>%
      group_by(sample_id) %>%
      dplyr::filter(!any(variant %in% EXCLUDE_SET[["BOTH"]]) & !(drug_short == "KAN" & sum(variant %in% EXCLUDE_SET[["KAN_EXT"]]) > 1)) %>%
      dplyr::filter(any(variant %in% EXCLUDE_SET[[paste0(drugName, "_EXT")]] & !het_strict)) %>%
      assignEpiGroup("eis") %>%
      dplyr::filter(variant %in% EXCLUDE_SET[[paste0(drugName, "_EXT")]] & !het_strict) %>%
      ungroup() %>%
      summarizeEpiGroup(c("Group", "variant", "phenotype"), "drug")
    epiTabs[[drugName]] = curEpiTab
  }
  curEpiTab = fullDataset %>%
    dplyr::filter(drug_short == "BDQ") %>%
    group_by(sample_id) %>%
    dplyr::filter(!any(gene %in% EXCLUDE_BDQ_GENES & Final <= 2)) %>%
    dplyr::filter(any(gene %in% BDQ_GENE & Final <= 2 & !het_strict))
  for (geneName in STRATIFY_BDQ_GENES) {
    subEpiTab = curEpiTab %>%
      assignEpiGroup(geneName) %>%
      dplyr::filter(gene %in% BDQ_GENE & Final <= 2 & !het_strict) %>%
      slice(1) %>%
      ungroup() %>%
      summarizeEpiGroup(c("Group", "phenotype"), c("gene", "drug"))
    epiTabs[[paste0("BDQ_", geneName)]] = subEpiTab
  }
  fullDataset %<>% group_by(sample_id, drug)
  list(fullDataset = fullDataset, epiTabs = epiTabs)
}

#' Compute compensatory mutation tables for INH
#' @noRd
computeCompensatoryStats = function(fullDataset) {
  fullDataset %<>% ungroup()
  compTabs = vector("list", 2) %>%
    magrittr::set_names(c("INH", "RIF"))
  miniTab1 = fullDataset %>%
    dplyr::filter(drug_short == "INH") %>%
    dplyr::filter(!(gene %in% c("ahpC", "katG") & Final >= 4)) %>%
    group_by(sample_id) %>%
    dplyr::filter(!any(is.na(phenotype))) %>%
    dplyr::filter(!any(gene == "ahpC" & effect %in% POOLED_EFFECTS[["LoF"]])) %>%
    ungroup() %>%
    mutate(MOI = (gene == "ahpC" & (!is.na(as.integer(position))) & as.integer(position) >= 2726101 & as.integer(position) <= 2726192)) %>%
    group_by(sample_id) %>%
    mutate(num_MOI = sum(MOI, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(relevant_MOI = ifelse(MOI & num_MOI == 1 & !het_strict, mutation, NA)) %>%
    group_by(sample_id) %>%
    dplyr::filter(any(!is.na(relevant_MOI))) %>%
    mutate(relevant_MOI = na.omit(unique(relevant_MOI)))
  allTabs = vector("list", 6) %>%
    magrittr::set_names(LETTERS[1:6])
  allTabs[["A"]] = miniTab1
  allTabs[["B"]] = miniTab1 %<>% dplyr::filter(!any(gene == "inhA" & Final <= 2))
  allTabs[["C"]] = miniTab1 %<>% dplyr::filter(!any(variant %in% paste0("katG_p.Ser315", c("Arg", "Asn", "Gly", "Ile", "Thr"))))
  allTabs[["D"]] = miniTab1 %<>% dplyr::filter(!any(gene == "katG" & (Final <= 2 | effect %in% POOLED_EFFECTS[["LoF"]])))
  allTabs[["E"]] = miniTab1 %<>% dplyr::filter(!any(gene == "katG" & Final == 3 & !(effect %in% SILENT_EFFECTS | effect == "upstream_gene_variant")))
  allTabs[["F"]] = miniTab1 %<>% dplyr::filter(!any(gene == "katG" & Final == 3 & !(effect %in% SILENT_EFFECTS)))
  compTabs[["INH"]] = allTabs
  fullDataset %<>% group_by(sample_id, drug)
  list(fullDataset = fullDataset, compTabs = compTabs)
}

#' Compute sensitivity and specificity statistics for a catalogue against a dataset
#' @param fullDataset A data frame of genotype-phenotype data as produced by loadSensSpecData.
#' @param safe Logical; if TRUE, applies the neutral algorithm from scratch before computing sens/spec.
#' @param skipEpistasis Logical; if TRUE, skips the epistasis grouping step.
#' @param sameRIF Logical; if TRUE, restricts sens/spec computation to isolates with RIF phenotype data.
#' @param minQ Minimum quality threshold to apply; NA skips quality filtering.
#' @param relaxed Logical; if TRUE, uses the relaxed MAF threshold for het classification.
#' @param skipCompensatory Logical; if TRUE, skips the compensatory mutation analysis.
#' @param saveIntermed Logical; if TRUE, writes intermediate result tables to disk.
#' @inheritParams applyCatalogue
#' @export
computeSensSpec = function(fullDataset,
                           catalogueFile,
                           version = CURR_VERSION,
                           safe = TRUE,
                           skipEpistasis = TRUE,
                           sameRIF = TRUE,
                           minQ = NA,
                           relaxed = FALSE,
                           skipCompensatory = TRUE,
                           saveIntermed = FALSE) {
  processedDatasets = imap(MAF_THRESHOLDS, function(threshold, curName) {
    print(paste("Processing option", curName))
    fullDataset %>%
      applyCatalogue(catalogueFile, minMAF = threshold, minQ = QUALITY_THRESHOLD_STRICT, LoF = TRUE, version = version)
  })
  fullDataset = processedDatasets[["regular"]] %>%
    mutate(het_relaxed = processedDatasets[["relaxed"]]$het) %>%
    mutate(het_strict  = processedDatasets[["strict" ]]$het)
  fullDataset = fullDataset %>%
    select(sample_id, drug, variant, phenotype, gene, mutation, effect, pos1, pos2, position, 
           het, het_relaxed, het_strict, Final, RRDR_NON_SILENT, LoF_candidate, lineage)
  ## Calculate the RIF resistance variable
  fullDataset = fullDataset %>%
    mutate(drug_short = SHORT_NAMES[match(drug, DRUG_LIST)]) %>%
    group_by(sample_id) %>%
    mutate(genoRIF          = (drug_short == "RIF" & ((!is.na(Final) & Final <= 2) | RRDR_NON_SILENT))) %>%
    mutate(RIF_geno_Relaxed = any(genoRIF & !het_relaxed), RIF_geno_Strict  = any(genoRIF & !het_strict))
  if (sameRIF) {
    fullDataset %<>%
      mutate(RIF_geno_Regular = RIF_geno_Relaxed)
  } else {
    fullDataset %<>%
      mutate(RIF_geno_Regular = any(genoRIF & !het))
  }
  fullDataset %<>% 
    ungroup()
  fullDataset %<>%
    mutate(Final_Relaxed = Final) %>%
    mutate_at("Final",         ~{ifelse(is.na(.)  & !het         & (RRDR_NON_SILENT | LoF_candidate), 2, .)}) %>%
    mutate_at("Final_Relaxed", ~{ifelse(is.na(.)  & !het_relaxed & (RRDR_NON_SILENT | LoF_candidate), 2, .)}) %>%
    mutate_at(c("Final", "Final_Relaxed"), ~{replace_na(., 3)}) %>%
    select(-pos1, -pos2, -LoF_candidate)
  fullDataset %<>%
    group_by(sample_id, drug)
  ## TODO: I still need to implement the epistasis rules for AMI, KAN and BDQ/CFZ in version 3!
  if (version == CURR_VERSION) {
    if (!skipEpistasis) {
      result      = computeEpistasisStats(fullDataset)
      fullDataset = result$fullDataset
      epiTabs     = result$epiTabs
    }
    if (!skipCompensatory) {
      result      = computeCompensatoryStats(fullDataset)
      fullDataset = result$fullDataset
      compTabs    = result$compTabs
    }
  }
  ## Assign the final "regular" and "relaxed" groups to each sample; NB: samples flagged for epistasis will be counted as not fitting the catalogue criteria
  ## The computation below adds MAX_GRADE to any het variant so that any group of interest (1, 2 or 3) can only get determined by relevant non-hets
  # if (version == CURR_VERSION) { ## NOTE: mark epistasis candidates as MAX_GRADE first, then take min
  #   fullDataset = fullDataset %>%
  #     mutate_at("het", ~{ifelse(excludeEpi_Regular, TRUE, .)}) %>%
  #     mutate_at("het_relaxed", ~{ifelse(excludeEpi_Relaxed, TRUE, .)})
  # }
  fullDataset = fullDataset %>%
    mutate(Group_Regular = min(Final + het * MAX_GRADE), Group_Relaxed = min(Final_Relaxed + het_relaxed * MAX_GRADE)) %>%
    mutate_at(c("Group_Regular", "Group_Relaxed"), ~{ifelse(. >= 3, . + 1, .)}) %>%
    mutate(extendedCandidate = (!is.na(effect) & effect != "upstream_gene_variant" & drug_short == "INH" & gene == "katG" & Final == 3)) %>%
    mutate_at("Group_Regular", ~{ifelse(any(variant %in% COMPENSATORY) & any((!het)         & extendedCandidate & . == 4), 3, .)}) %>%
    mutate_at("Group_Relaxed", ~{ifelse(any(variant %in% COMPENSATORY) & any((!het_relaxed) & extendedCandidate & . == 4), 3, .)}) %>%
    select(-extendedCandidate)
    # mutate(extendedCandidate = (effect %in% ADDITIONAL_EFFECTS & (gene %in% unlist(ADDITIONAL_GENES[-1]) | (drug_short == "DLM" & gene %in% EXTENDED_ADD_GENES[["DLM"]])))) %>%
    # mutate_at("Group_Regular", ~{ifelse(RIF_geno_Regular & any(!het         & Final == 3 & extendedCandidate & . == 4), 3, .)}) %>%
    # mutate_at("Group_Relaxed", ~{ifelse(RIF_geno_Relaxed & any(!het_relaxed & Final == 3 & extendedCandidate & . == 4), 3, .)}) %>%
    # select(-extendedCandidate)
  if (safe) { 
    stopifnot(all(testConsistent(fullDataset, groupingVars = c("sample_id", "drug"), consistentVars = c("Group_Regular", "Group_Relaxed"))[[1]])) 
    if (useLineageData) {
      stopifnot(all(testConsistent(fullDataset, groupingVars = c("sample_id"), consistentVars = c("lineage"))[[1]])) 
    }
  }
  fullDataset = fullDataset %>%
    slice(1) %>%
    ungroup()
  ## Calculate the size of groups 1, 2, 3, 1 + 2, and 1 + 2 + 3, with both regular and relaxed thresholds and each RIF_geno status, for each drug
  finalTabs = vector("list", 3 * 2 * 5 * 8) %>%
    magrittr::set_names(outer(paste0("Lineage", c(1:6, "Other", "Any")), outer(c("R", "S", "ALL"), outer(c("Regular", "Relaxed"), STAGES, 
                    function(x, y) {paste0(x, "_", y)}), function(z, w) {paste0("RIF", z, "_", w)}), function(a, b) {paste0(a, "_", b)}))
  for (relax in c(FALSE, TRUE)) {
    if (relax) {
      curDataset = fullDataset %>% 
        mutate(relevantGroup = Group_Relaxed, relevantGeno = RIF_geno_Relaxed)
    } else {
      curDataset = fullDataset %>% 
        mutate(relevantGroup = Group_Regular, relevantGeno = RIF_geno_Regular)
    }
    for (geno in c("R", "S", "ALL")) {
      if (geno != "ALL") {
        curSubset = curDataset %>%
          dplyr::filter(relevantGeno == ifelse(geno == "R", TRUE, FALSE))
      } else {
        curSubset = curDataset
      }
      for (Lineage in c("Any", as.character(1:6), "Other")) {
        if (Lineage == "Any") {
          curSuperTab = curSubset
        } else {
          curSuperTab = curSubset %>%
            dplyr::filter(lineage == Lineage)
        }
        for (stage in STAGES) {
          if (stage <= 3) {
            curTab = curSuperTab %>%
              mutate(selected = (relevantGroup == stage & (stage <= 2 | (stage == 3 & drug %in% GROUP_3_DRUGS))))
          } else {
            if (stage == 12) {
              curTab = curSuperTab %>%
                mutate(selected = (relevantGroup <= 2))
            } else { ## stage = 123
              curTab = curSuperTab %>%
                mutate(selected = (relevantGroup <= 2 | (relevantGroup == 3 & drug %in% GROUP_3_DRUGS)))
            }
          }
          curTab %<>% 
            mutate_at("selected", ~{convertToLogical(.)}) %>%
            ungroup()
          if (saveIntermed) {
            write_csv(curTab, file = paste0("SensSpec_Intermediate_RIF", geno, "_", ifelse(relax, "Relaxed", "Regular"), "_Group", stage, "_Lineage", Lineage, ".csv"))
          }
          curTab %<>%
            group_by(selected, phenotype, drug) %>%
            summarise(N = n_distinct(sample_id), .groups = "drop")
          finalTabs[[paste0("Lineage", Lineage, "_", "RIF", geno, "_", ifelse(relax, "Relaxed", "Regular"), "_", stage)]] = curTab
        }
      }
    }
  }
  output = list(finalTabs = finalTabs)
  if (!skipEpistasis) {
    output = c(output, list(epiTabs = epiTabs))  
  }
  if (!skipCompensatory) {
    output = c(output, list(compTabs = compTabs))
  }
  postprocessTabs(output, version = version, relaxed = relaxed, sameRIF = sameRIF, minQ = minQ)
}

#' Post-process sensitivity/specificity result tables into a final summary format
#' @noRd
postprocessTabs = function(List, version, relaxed, sameRIF, minQ) {
  fullTab = tibble()
  output = c()
  if ("finalTabs" %in% names(List)) {
    longTabs = List[["finalTabs"]]
    for (curName in names(longTabs)) {
      curTab = longTabs[[curName]]
      if (nrow(curTab) == 0) next
      curTab %<>%
        group_by(drug) %>%
        summarise(TP = max(N * as.integer(selected & phenotype == "R")), TN = max(N * as.integer(!selected & phenotype == "S")),
                  FP = max(N * as.integer(selected & phenotype == "S")), FN = max(N * as.integer(!selected & phenotype == "R")),
                  .groups = "drop") %>%
        mutate(group = curName, .before = 1)
      fullTab %<>%
        bind_rows(curTab)
    }
    fullTab %<>%
      mutate(PPV  = map2(TP, TP + FP, safeBinomTest), NPV  = map2(TN, TN + FN, safeBinomTest),
             Sens = map2(TP, TP + FN, safeBinomTest), Spec = map2(TN, FP + TN, safeBinomTest), propR = map2(TP + FN, TP + FN + FP + TN, safeBinomTest)) %>%
      expandBinomCIs()
    if (!is.na(minQ)) {
      fullTab %<>%
        dplyr::filter(str_detect(group, ifelse(relaxed, "Relaxed", "Regular")))
    }
    output = c(output, list(fullTab = fullTab))
  }
  if ("epiTabs" %in% names(List)) {
    shortTabs = List[["epiTabs"]]
    miniTab   = tibble()
    for (curName in names(shortTabs)) {
      curTab = shortTabs[[curName]]
      if ("variant" %in% colnames(curTab)) {
        curTab %<>% rename(criterion = variant)
      } else {
        curTab %<>% rename(criterion = gene) 
      }
      curTab %<>%
        group_by(criterion) %>%
        summarise(LOF_R   = max(N * as.integer(Group == "A" & phenotype == "R")), LOF_S   = max(N * as.integer(Group == "A" & phenotype == "S")),
                  nomut_R = max(N * as.integer(Group == "B" & phenotype == "R")), nomut_S = max(N * as.integer(Group == "B" & phenotype == "S")),
                  .groups = "drop") %>%
        mutate(group = curName, .before = 1)
      miniTab %<>%
        bind_rows(curTab)
    }
    miniTab %<>%
      mutate(PPV_LOF = map2(LOF_R, LOF_R + LOF_S, safeBinomTest), PPV_nomut  = map2(nomut_R, nomut_R + nomut_S, safeBinomTest)) %>%
      expandBinomCIs()
    output = c(output, list(epiTab = miniTab))
  }
  if ("compTabs" %in% names(List)) {
    extraTabs = List[["compTabs"]]
    miniTab   = tibble()
    for (curDrug in names(extraTabs)) {
      fullTab = extraTabs[[curDrug]]
      if (curDrug == "INH") {
        for (curName in names(fullTab)) {
          curTab = fullTab[[curName]] %>%
            mutate(mut_R = (phenotype == "R"), mut_S = (phenotype == "S"), 
                   mut_R_genoRIF = (phenotype == "R" & RIF_geno_Strict), mut_S_genoRIF = (phenotype == "S" & RIF_geno_Strict)) %>%
            slice(1) %>%
            ungroup() %>%
            group_by(relevant_MOI) %>%
            summarise(across(c(mut_R, mut_S, mut_R_genoRIF, mut_S_genoRIF), sum), .groups = "drop") %>%
            mutate(drug = curDrug, group = curName, .before = 1)
          miniTab %<>%
            bind_rows(curTab)
        }
      }
    }
    miniTab %<>%
      mutate(PPV = map2(mut_R, mut_R + mut_S, safeBinomTest)) %>%
      expandBinomCIs()
    output = c(output, list(compTab = miniTab))
  }
  output
}

#' Write sensitivity/specificity results to CSV files
#' @noRd
writeSensSpecOutputs = function(result, version, relaxed, sameRIF, minQ) {
  if ("fullTab" %in% names(result)) {
    write_csv(result$fullTab, paste0("SensSpec_Leonid_Version", version, "_", ifelse(relaxed, "Relaxed_", ""),
                                     ifelse(sameRIF, "sameRIF", ""), "_", ifelse(!is.na(minQ), paste0("Q", minQ), ""), "_", Sys.Date(), ".csv"))
  }
  if ("epiTab" %in% names(result)) {
    write_csv(result$epiTab, paste0("Epistasis_Leonid_Version", version, "_", ifelse(relaxed, "Relaxed_", ""), Sys.Date(), ".csv"))
  }
  if ("compTab" %in% names(result)) {
    write_csv(result$compTab, paste0("Compensatory_Leonid_Version", version, "_", ifelse(relaxed, "Relaxed_", ""), Sys.Date(), ".csv"))
  }
}

#' Run the full sensitivity/specificity computation pipeline
#' @noRd
runComputeSensSpec = function(version = CURR_VERSION,
                              safe = TRUE,
                              skipEpistasis = TRUE,
                              sameRIF = TRUE,
                              minQ = NA,
                              relaxed = FALSE,
                              skipCompensatory = TRUE,
                              useOrphanData = TRUE,
                              useLineageData = TRUE,
                              useSublineageData = FALSE,
                              saveIntermed = FALSE,
                              catalogueFile = NULL,
                              DATA_DIRECTORY = "SOLO Algorithm Input files/DATABASE EXTRACTION files/") {
  if (is.null(catalogueFile)) {
    if (version == CURR_VERSION) {
      catalogueFile = paste0("Results/", EXTRACTION_ID, "/List_of_graded_variants_", format(Sys.Date(), "%d%b%Y"), ".csv")
    } else {
      catalogueFile = "who_catalogue_prev_version.xlsx"
    }
  }
  fullDataset = loadSensSpecData(DATA_DIRECTORY, version, useOrphanData, useLineageData, useSublineageData)
  result = computeSensSpec(fullDataset, catalogueFile, version, safe, skipEpistasis, sameRIF, minQ, relaxed, skipCompensatory, saveIntermed)
  writeSensSpecOutputs(result, version, relaxed, sameRIF, minQ)
  result
}

#' Merge version 2 and version 3 sensitivity/specificity output files
#' @noRd
mergeOutputs = function(file2 = "SensSpec_Leonid_Version2_sameRIF__2026-02-25.csv",
                        file3 = "SensSpec_Leonid_Version3_sameRIF__2026-02-25.csv") {
  tab2 = read_csv(file2, show_col_types = FALSE)
  tab3 = read_csv(file3, show_col_types = FALSE)
  mergedTab = tab2 %>%
    inner_join(tab3, by = c("group", "drug"), suffix = c("_v2", "_v3"))
  cnames = colnames(tab2)
  cnames2 = paste0(cnames[-(1:2)], "_v2")
  cnames3 = paste0(cnames[-(1:2)], "_v3")
  cnames = as.vector(rbind(cnames2, cnames3))
  mergedTab = mergedTab %>%
    select(all_of(c("group", "drug", cnames)))
  write_csv(mergedTab, paste0("SensSpec_Leonid_Version2_vs_Version3_", Sys.Date(), ".csv"))
  mergedTab
}
