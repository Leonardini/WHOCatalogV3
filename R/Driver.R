#' Preprocess raw genotype data by renaming columns and applying MAF and quality filters
#' @noRd
preprocessGenotypes = function(allGenotypes, minMAF = MAF_THRESHOLD_REGULAR, lowMAFHet = TRUE, minQ = QUALITY_THRESHOLD_STRICT, lowQHet = TRUE) {
  ## Rename three columns for convenience and remove an unused column
  allGenotypes %>%
    rename(gene = 'resolved_symbol', mutation = 'variant_category', effect = 'predicted_effect') %>%
    select(-neutral) %>%
    ## Replace missing mutations or effects by "missing", then create a variant out of gene and mutation
    mutate_at(c("mutation", "effect"), ~{replace_na(., "missing")}) %>%
    mutate(variant = ifelse(mutation == "missing", "missing", paste(gene, mutation, sep = "_"))) %>%
    ## Mark any genotype entry with a low MAF or a missing variant (i.e. a sequencing defect in the gene) as a het
    applyThresholds(minMAF = minMAF, lowMAFHet = lowMAFHet, minQ = minQ, lowQHet = lowQHet) %>%
    extractPositions(colName = "mutation", maxNumber = 2L) %>%
    ## Marking non-silent variants in the RRDR_GENE with at least one of the two AAs in the RRDR (originally this was in Neutral.R)
    computeRRDRInfo()
}

#' Identify samples that fail quality control based on specific variant-drug combinations
#' @noRd
identifySamplesToExclude = function(mainDataset) {
  ## Quality control based on a sanity check: S isolates with one of specific variant-drug combinations fail QC
  ## For bookkeeping purposes we keep additional columns even though only the sample_ids are used for exclusion
  mainDataset %>%
    inner_join(BAD_VAR_DRUG_PAIRS, by = c("drug", "variant")) %>%
    dplyr::filter(phenotype == "S") %>%
    mutate(neutral = FALSE) %>%
    distinct(sample_id, drug, variant, gene, phenotype, het, effect, tier, neutral, category_phenotype)
}

#' Annotate all datasets with neutral variants, exclusions, and denominators
#' @noRd
annotateDatasets = function(fullDataset, samplesToExclude, allNeutrals) {
  for (name in names(fullDataset)) {
    if (name %in% c("MAIN", "WHO")) {
      print(paste("Preprocessing the", name, "dataset"))
      curSet = fullDataset[[name]]
      ## Compute the per-drug denominators (number of R/S isolates screened for each drug); account for exclusion
      allDenominators = curSet %>%
        dplyr::filter(!(sample_id %in% samplesToExclude$sample_id)) %>%
        group_by(drug) %>%
        distinct(sample_id, .keep_all = TRUE) %>%
        mutate(RDen = sum(phenotype == "R"), SDen = sum(phenotype == "S")) %>%
        slice(1) %>%
        ungroup() %>%
        select(drug, RDen, SDen)
      curSet %<>%
        inner_join(allDenominators, by = "drug")
      ## Adjust for the numbers of R/S isolates screened for each drug but excluded due to QC just for the bad pairs
      extraDenominators = curSet %>%
        dplyr::filter(sample_id %in% samplesToExclude$sample_id) %>%
        group_by(drug) %>%
        distinct(sample_id, .keep_all = TRUE) %>%
        mutate(RDen = sum(phenotype == "R"), SDen = sum(phenotype == "S")) %>%
        slice(1) %>%
        ungroup() %>%
        select(-variant) %>%
        inner_join(BAD_VAR_DRUG_PAIRS, by = "drug") %>%
        select(drug, variant, RDen, SDen)
      curSet %<>%
        left_join(extraDenominators, by = c("drug", "variant")) %>%
        adjustDuplicateColumns(warn = FALSE, add = TRUE)
      ## Mark the neutral mutations and convert the neutral variable to a logical one
      curSet %<>%
        left_join(allNeutrals, by = c("drug", "variant")) %>%
        mutate_at(c(paste0("set", LETTERS[1:5]), "lit_mutation", "prev_version", "neutral"), convertToLogical)
      ## And finally, replace into the full dataset
      fullDataset[[name]] = curSet
    }
  }
  fullDataset
}

#' Compute per-variant frequency and solo counts for QC-excluded samples
#' @noRd
computeExcludedCounts = function(curSet, samplesToExclude, datasetName) {
  ## Exclude the samples failing QC, after saving the rows corresponding to the problematic mutations first
  ## The samples that were among the originally excluded ones are accounted for as well for completeness
  curExcludedEntries = curSet %>%
    dplyr::filter(sample_id %in% samplesToExclude$sample_id) %>%
    bind_rows(samplesToExclude %>% dplyr::filter(datasetName == "MAIN" | (category_phenotype == datasetName))) %>%
    distinct(sample_id, drug, variant, gene, phenotype, het, effect, tier, neutral)
  ## Count frequency and solo frequency among the excluded samples, in the entries with "QC-relevant" mutations
  curExtraCounts = curExcludedEntries %>%
    dplyr::filter(!het) %>%
    inner_join(BAD_VAR_DRUG_PAIRS, by = c("drug", "variant")) %>%
    group_by(drug, variant) %>%
    mutate(present_R = sum(phenotype == "R"), present_S = sum(phenotype == "S")) %>%
    slice(1) %>%
    ungroup() %>%
    select(drug, variant, present_R, present_S)
  curExtraOutputs = curExcludedEntries %>%
    prepMask(Silent = TRUE, Tier2 = TRUE, Neutral = TRUE, Pool = NA, SOnly = TRUE) %>%
    runSOLOPipeline(maxIter = 1L, stage = NULL, removeSOnly = FALSE)
  curExtraSolos = curExtraOutputs[[1]] %>%
    inner_join(BAD_VAR_DRUG_PAIRS, by = c("drug", "variant")) %>%
    select(drug, variant, Rcnt, Scnt)
  list(extraCounts = curExtraCounts, extraSolos = curExtraSolos)
}

#' Apply manual grading overrides to the final catalogue
#' @noRd
applyManualChecks = function(finalCatalog, manualCheckResults) {
  ## Replace the grading and the supplementary considerations for the variants that were manually checked; this is a side effect by design
  finalCatalog %>%
    full_join(manualCheckResults, by = c("drug", "variant"), suffix = c("_auto", "_manual")) %>%
    adjustDuplicateColumns(warn = FALSE, add = FALSE, suffixes = c("_manual", "_auto"))
}

#' Compute per-variant frequency and SOLO statistics, adjusted for QC-excluded samples
#' @noRd
computeVariantStats = function(curSet, curExtraCounts, curExtraSolos) {
  ## Calculate the frequency of each variant in R isolates (present_R), in S isolates (present_S), and overall
  ## (present), and corresponding variables for the complement (absent), as well as SOLO counts per phenotype.
  ## Note that hets are removed at the beginning!
  ## Also adjust the variant and solo counts for samples excluded due to QC by adding the extra counts to each.
  curStats = curSet %>%
    dplyr::filter(!het & variant != "missing") %>% ## TODO: REVISIT THE LOGIC OF MISSING VARIANTS IN SENS-SPEC!
    rename(datasets = category_phenotype) %>%
    group_by(drug, variant) %>%
    mutate(present = n(), present_R = sum(phenotype == "R"), present_S = present - present_R) %>%
    left_join(curExtraCounts, by = c("drug", "variant")) %>%
    adjustDuplicateColumns(warn = FALSE, add = TRUE) %>%
    mutate(absent_R = RDen - present_R, absent_S = SDen - present_S) %>%
    left_join(curExtraSolos, by = c("drug", "variant")) %>%
    adjustDuplicateColumns(warn = FALSE, add = TRUE) %>%
    ungroup() %>%
    distinct(drug, variant, .keep_all = TRUE)
  ## Finally, complete the extraction and mark variants to FDR-correct for in SOLO and ALL modes, respectively.
  curStats %>%
    mutate(correctAll = (!(lit_mutation | prev_version | is.na(effect) | effect %in% SILENT_EFFECTS))) %>%
    rename(SOLO_R = Rcnt, SOLO_S = Scnt) %>%
    mutate(SOLO_SorR = SOLO_R + SOLO_S, correctSOLO = correctAll & !is.na(SOLO_SorR) & SOLO_SorR > 0) %>%
    select(drug, variant, stage, tier, neutral, datasets, effect, position, stage, lit_mutation, prev_version,
           starts_with("present"), starts_with("absent"), starts_with("SOLO"), starts_with("correct"), starts_with("pos"), starts_with("set"))
}

#' Adjust variant statistics after removing specified mutations and their associated isolates
#' @noRd
adjustStatsForRemovedMutations = function(curStats, curDataset, mutationsToRemove, curStage, soloIsolates) {
  extraIsolates = curDataset %>%
    inner_join(mutationsToRemove, by = c("drug", "variant")) %>%
    pull(sample_id) %>%
    unique()
  print(paste("Removing", length(extraIsolates), "isolates at stage", curStage))
  extraDenominators = curDataset %>% ## TODO: Make sure that we only need to subtract non-het non-missing counts!
    dplyr::filter(stage == curStage & sample_id %in% extraIsolates & !het & variant != "missing") %>%
    group_by(drug, variant) %>%
    mutate(present = n(), present_R = sum(phenotype == "R"), present_S = present - present_R) %>%
    mutate_at(paste0('present', c('', '_R', '_S')), ~{multiply_by(., -1)}) %>%
    slice(1) %>%
    ungroup() %>%
    select(drug, variant, starts_with('present'))
  # Calculate SOLO counts to subtract (NEW) - suggested by Claude!
  extraSoloSubtract = tibble(drug = character(), variant = character(), SOLO_R = integer(), SOLO_S = integer())
  if (length(extraIsolates) > 0) {
    extraPhenotypes = curDataset %>%
      dplyr::filter(sample_id %in% extraIsolates & stage == curStage & variant != "missing") %>%
      distinct(drug, sample_id, phenotype)
    extraSoloSubtract = soloIsolates %>%
      dplyr::filter(sample_id %in% extraIsolates) %>%
      inner_join(extraPhenotypes, by = c("drug", "sample_id")) %>%
      group_by(drug, variant) %>%
      summarise(SOLO_R = -sum(phenotype == "R"), SOLO_S = -sum(phenotype == "S"), .groups = "drop")
  }
  # Fix suggested by Claude
  curStats %>%
    mutate(RDen = present_R + absent_R, SDen = present_S + absent_S) %>%
    left_join(extraDenominators, by = c("drug", "variant")) %>%
    adjustDuplicateColumns(warn = FALSE, add = TRUE) %>%
    left_join(extraSoloSubtract, by = c("drug", "variant")) %>%  # NEW
    adjustDuplicateColumns(warn = FALSE, add = TRUE) %>%          # NEW
    mutate(absent_R = RDen - present_R, absent_S = SDen - present_S) %>%
    mutate(SOLO_SorR = SOLO_R + SOLO_S) %>%  # NEW: recalculate
    select(-RDen, -SDen) %>%
    distinct(drug, variant, .keep_all = TRUE)
}

#' Execute the complete analysis required to create version 2 or 3 of the WHO catalogue
#' @param correct_all Logical; if TRUE, FDR correction uses all variants; if FALSE, only variants appearing as a SOLO in at least one isolate.
#' @param LoF Logical; if TRUE, pooled LoF mutations are included as hypotheses for FDR correction.
#' @param safe Logical; if TRUE, performs a full from-scratch conversion of variants across versions (slow).
#' @param minMAF Minimum minor allele frequency; variants below this threshold are treated as heterozygous.
#' @param lowMAFHet Logical; if TRUE, low-MAF variants become hets; if FALSE, they are filtered out.
#' @param minQ Minimum quality score; variants below this threshold are treated as heterozygous.
#' @param lowQHet Logical; if TRUE, low-quality variants become hets; if FALSE, they are filtered out.
#' @param EXTRACTION_ID Identifier string for the database extraction (used to locate input files); defaults to the package constant EXTRACTION_ID.
#' @param listIsolates Logical; if TRUE, records all isolates classified as solo at each iteration.
#' @param augmentWithOrphansAndLineage Logical; if TRUE, add a breakdown by sub-lineage and the number of orphans.
#' @param OUTPUT_DIRECTORY Directory for writing output files; defaults to Results/<EXTRACTION_ID>.
#' @param DATA_DIRECTORY Directory containing the database extraction files.
#' @param NON_DATABASE_DIRECTORY Directory containing non-database reference files (defaults to inst/extdata).
#' @export
mainDriver = function(correct_all = TRUE,
                      LoF = TRUE, 
                      safe = TRUE,
                      minMAF = MAF_THRESHOLD_REGULAR, 
                      lowMAFHet = TRUE, 
                      minQ = QUALITY_THRESHOLD_STRICT, 
                      lowQHet = TRUE, 
                      EXTRACTION_ID = EXTRACTION_ID,
                      listIsolates = TRUE,
                      augmentWithOrphansAndLineage = TRUE,
                      OUTPUT_DIRECTORY = paste0("Results/", EXTRACTION_ID),
                      DATA_DIRECTORY         = "SOLO Algorithm Input files/DATABASE EXTRACTION files/",
                      NON_DATABASE_DIRECTORY = system.file("extdata", package = "SOLOport")) {
  DATA_DIRECTORY %<>%
    normalizePath()
  NON_DATABASE_DIRECTORY %<>%
    normalizePath()
  dir.create(OUTPUT_DIRECTORY, recursive = TRUE)
  OUTPUT_DIRECTORY %<>%
    normalizePath()
  ## Initial preprocessing of the genotypes: extract all the genotypes, recording the drug and the tier
  print("Extracting the genotypes")
  allGenotypes  = extractData(inDir = paste(str_remove(DATA_DIRECTORY, "/$"), str_remove(EXTRACTION_ID, "/$"), "full_genotypes/", sep = "/"),
                              drugList = DRUG_LIST, geno = TRUE)
  ## Make sure there are no missing or empty resolved symbols
  stopifnot(!any(is.na(allGenotypes$resolved_symbol)))
  allGenotypes = preprocessGenotypes(allGenotypes, minMAF = minMAF, minQ = minQ)
  runGenotypeConsistencyTests(allGenotypes)
  ## Initial preprocessing of the phenotypes: extract all the phenotypes, recording the drug
  print("Extracting the phenotypes")
  allPhenotypes = extractData(inDir = paste(str_remove(DATA_DIRECTORY, "/$"), str_remove(EXTRACTION_ID, "/$"), "phenotypes/", sep="/"),
                              drugList = DRUG_LIST, geno = FALSE)
  ## Rename one column for convenience and remove an unused column
  allPhenotypes = allPhenotypes %>% 
    rename(category_phenotype = 'phenotypic_category') %>%
    select(-"box")
  runPhenotypeConsistencyTests(allPhenotypes)
  ## Merge the genotypes and the phenotypes after separating the latter by phenotypic category group
  fullDataset = mergeGenoPheno(allGenotypes, allPhenotypes, phenoGroups = PHENO_GROUPS)
  samplesToExclude = identifySamplesToExclude(fullDataset[["MAIN"]])
  ## Record the sample_ids that need to be excluded
  write_csv(samplesToExclude, file.path(OUTPUT_DIRECTORY, "excluded_after_qc.csv"))
  ## Set the data with WHO phenotype and add it to the dataset; it will be used to determine neutral variants
  WHOData = fullDataset[["MAIN"]] %>%
    dplyr::filter(category_phenotype == "WHO")
  fullDataset[["WHO"]] = WHOData
  ## Replace the phenotype in the MAIN dataset with "ALL"
  fullDataset[["MAIN"]] %<>%
    mutate(category_phenotype = "ALL")
  ## (Re)run the neutral algorithm to create the files containing neutral mutations (note that this is a side effect)
  if (safe) {
    litTab = read_csv(file.path(str_remove(NON_DATABASE_DIRECTORY, "/$"), "literature_neutrals.csv"),
                      show_col_types = FALSE, guess_max = LARGE_NUMBER)
    WHODataMarked = WHOData %>%
      dplyr::filter(!(sample_id %in% samplesToExclude$sample_id)) %>%
      neutralAlgorithm(litTab = litTab, NON_DATABASE_DIRECTORY = str_remove(NON_DATABASE_DIRECTORY, "/$"))
    writeNeutralOutputs(WHODataMarked, safe = safe, outDir = OUTPUT_DIRECTORY)
  }
  ## Extract the prepared list of all neutral mutations
  ## TODO: Consider doing this with an in-memory variable
  allNeutrals = read_csv(file.path(OUTPUT_DIRECTORY, "neutral_mutations_WHO_F.csv"), guess_max = LARGE_NUMBER, show_col_types = FALSE) %>%
    mutate(neutral = TRUE)
  ## Now identify the neutral variants and mark these, as well as other relevant ones, in the entire dataset
  fullDataset = annotateDatasets(fullDataset, samplesToExclude, allNeutrals)
  ## Save the clean version of the input data for downstream analysis of the catalogue's performance on itself 
  write_csv(fullDataset[["MAIN"]], file.path(OUTPUT_DIRECTORY, "CompleteDataset.csv"))
  write_csv(fullDataset[["WHO"]], file.path(OUTPUT_DIRECTORY, "CompleteDatasetWHO.csv"))
  stageStats = vector("list", nrow(PHENO_GROUPS) * (length(POOLED_EFFECTS) + 1))
  ## Now process the entire dataset in the specified order; additionally perform pooling
  for (name in c("WHO", "MAIN")) { ## SEPT 8, 2025: Removing the CC and ATU options as they were only used in v2
    ## Placeholder for the primary statistics (and possibly the solo isolates)
    fullStats = tibble()
    fullIsolates = tibble()
    ## Sanity check: ensures that there is only one non-missing variant of each name per drug-sample combination
    stopifnot(anyDuplicated(fullDataset[[name]] %>% dplyr::filter(variant != "missing") %>% select(drug, sample_id, variant)) == 0)
    for (pool in c(NA, names(POOLED_EFFECTS))) {
      curIsolates = tibble()
      POOLED = !is.na(pool)
      print(paste("Processing the combination of", name, "and", ifelse(POOLED, pool, "no"), "pooling"))
      ## Extract the dataset corresponding to the name
      curSet = fullDataset[[name]]
      excluded     = computeExcludedCounts(curSet, samplesToExclude, name)
      curExtraCounts = excluded$extraCounts
      curExtraSolos  = excluded$extraSolos
      curSet %<>%
        dplyr::filter(!(sample_id %in% samplesToExclude$sample_id))
      ## Add the pooling if necessary; otherwise, process stage 0 which yields SOLO counts for neutral variants
      if (POOLED) {
        ## Update the variants being pooled accordingly; compress identical variants per drug-sample-gene combo
        ## Update the rows of the variants being pooled with missing data in the variables that are undefined
        ## Create an additional pooled variant for all the RRDR_NON_SILENT mutations
        ## Consider all the variants in the pool to have a LoF effect
        curSet %<>%
          mutate_at("variant", ~{ ifelse(effect %in% POOLED_EFFECTS[[pool]] & !het, paste0(gene, "_", pool), .) }) %>%
          mutate_at("variant", ~{ ifelse(RRDR_NON_SILENT & !het, paste0("RRDR", "_", pool), .)}) %>%
          mutate_at(c("max(af)", "position", "max(quality)"), ~{ ifelse(str_ends(variant, pool), NA, .) }) %>%
          mutate_at("neutral", ~{ ifelse(str_ends(variant, pool), FALSE, .) }) %>%
          mutate_at("effect",  ~{ ifelse(str_ends(variant, pool), "LoF", .) }) %>%
          distinct(drug, sample_id, gene, variant, .keep_all = TRUE)
      } 
      curSet = markStages(curSet)
      for (stage in 1:ifelse(POOLED, 2, 3)) { ## NOTE: stage 0 is removed as it only counts the neutral variants
        print(paste("Processing stage", stage))
        if (stage == 1) {
          ## Prepare stage 1: mask all neutral, silent and tier 2 variants
          relevantSet = prepMask(curSet, Silent = TRUE,  Tier2 = TRUE,  Neutral = TRUE,  Pool = pool, SOnly = TRUE)
        } else if (stage == 2) {
          ## Prepare for stage 2: mask covers the neutral or silent effect variants, but not those in tier 2
          relevantSet = prepMask(curSet, Silent = TRUE,  Tier2 = FALSE, Neutral = TRUE,  Pool = pool, SOnly = TRUE)
        } else { ## stage = 3
          ## Prepare for stage 3: mask the neutral variants or those in tier 2, but unmask silent tier 1 variants
          relevantSet = prepMask(curSet, Silent = FALSE, Tier2 = TRUE,  Neutral = TRUE,  Pool = NA,   SOnly = TRUE)
        }
        ## Run the SOLO pipeline for 1 step of the algorithm; record stage
        relevantOutputs = runSOLOPipeline(relevantSet, maxIter = 1L, stage = stage, removeSOnly = TRUE, 
                                          listIsolates = listIsolates)
        if (listIsolates) {
          relevantIsolates = relevantOutputs[[2]]
          write_csv(relevantIsolates, file = file.path(OUTPUT_DIRECTORY, paste0('SOLOsAtStage', stage, '.csv')))
        }
        relevantOutputs = relevantOutputs[[1]]
        relevantOutputs = relevantOutputs %>%
          dplyr::filter(variant != "missing")
        if (stage == 1) { ## Save the first stage results for future reference
          firstOutputs = relevantOutputs
        } else { ## In stages 2 and 3, remove the results of any pairs that were analysed in the first stage (i.e. tier 1 variants)
          relevantOutputs %<>%
            anti_join(firstOutputs, by = c("drug", "variant"))
        }
        ## Remove the classification (we will only use SOLO counts), and merge the original data with the new SOLO counts
        curSet %<>%
          left_join(relevantOutputs %>% select(drug, variant, Rcnt, Scnt), by = c("drug", "variant")) %>%
          adjustDuplicateColumns()
      }
      curStats = computeVariantStats(curSet, curExtraCounts, curExtraSolos)
      ## Aggregate the primary statistics with the ones computed previously, but only add new variant-drug pairs
      if (nrow(fullStats) > 0) {
        curStats %<>%
          anti_join(fullStats, by = c("drug", "variant"))
      }
      stopifnot(!POOLED || all(str_ends(curStats$variant, pool)))
      if (POOLED) {
        curStats %<>%
          mutate(correctAll = LoF, correctSOLO = LoF)
      }
      fullStats %<>%
        bind_rows(curStats)
      ## Save the resulting dataset and the statistics (corresponding to stage 1)
      useName = paste0(name, "_", ifelse(!POOLED, "unpooled", pool))
      fullDataset[[useName]] = curSet
      stageStats[[useName]]  = curStats
    }
    ## Save the basic mutation-wise statistics into a file
    curFilename = paste0("Stats_", name, ifelse(LoF, "_withLoFs", ""), "_Stage1.csv")
    write_csv(fullStats, file.path(OUTPUT_DIRECTORY, paste0("Basic_", curFilename)))
    ## Now compute the derived statistics and save them into a file
    fullStats %<>%
      computeCatalogueStats(correct_all = correct_all)
    write_csv(fullStats, file.path(OUTPUT_DIRECTORY, curFilename))
  }
  print("Computing the final grades of all mutations")
  stopifnot(length(POOLED_EFFECTS) == 1) ## Ensure that there is exactly one pool, as otherwise the simplified logic below breaks!
  ## CHANGE ON SEPT 11, 2025: the grading is now stage-wise and samples with grade 1/2 mutations are ignored in subsequent stages!
  gradedCatalog = NULL
  for (curStage in 1:3) {
    print(paste("Currently processing stage", curStage))
    if (curStage > 1) {
      mutationsToRemove = gradedCatalog %>%
        dplyr::filter(Final %in% GRADES[1:2] & stage < curStage) %>%
        select(drug, variant)
      for (name in setdiff(names(fullDataset), unlist(PHENO_GROUPS))) {
        ## Count the number of instances of each mutation being lost
        curDataset = fullDataset[[name]]
        soloFile = file.path(OUTPUT_DIRECTORY, paste0('SOLOsAtStage', curStage, '.csv'))
        stopifnot(file.exists(soloFile))
        soloIsolates = read_csv(soloFile, show_col_types = FALSE, guess_max = LARGE_NUMBER)
        curStats = stageStats[[name]]
        adjustedStats = adjustStatsForRemovedMutations(curStats, curDataset, mutationsToRemove, curStage, soloIsolates)
        if (str_detect(name, "unpooled")) {
          initStats = adjustedStats
        } else {
          adjustedStats %<>%
            anti_join(initStats, by = c("drug", "variant")) %>%
            mutate(correctAll = LoF, correctSOLO = LoF)
          stopifnot(all(str_ends(adjustedStats$variant, pool)))
          initStats %<>%
            rbind(adjustedStats)
          if (str_ends(name, tail(names(POOLED_EFFECTS), 1))) {
            ## Compute the adjusted derived statistics and save them into a file
            initStats %<>%
              computeCatalogueStats(correct_all = correct_all)
            redName = str_split_fixed(name, "_", 2)[,1]
            curFilename = paste0("Stats_", redName, ifelse(LoF, "_withLoFs", ""), "_Stage", curStage, ".csv")
            write_csv(initStats, file.path(OUTPUT_DIRECTORY, curFilename))
          }
        }
      }
    }
    gradedCatalog = gradeMutations(LoF = LoF, NON_DATABASE_DIRECTORY = str_remove(NON_DATABASE_DIRECTORY, "/$"), stage = curStage, outDir = OUTPUT_DIRECTORY)
  }
  catalogFiles = file.path(OUTPUT_DIRECTORY, paste0("Final_graded_algorithm_catalogue_", Sys.Date(), "_withLoFs_Stage", 1:3, ".csv"))
  finalCatalog = tibble()
  for (curStage in 1:3) {
    curSubCatalog = read_csv(catalogFiles[curStage]) %>%
      dplyr::filter(stage == curStage | (curStage == 3 & (stage == 4 | is.na(stage)) & !is.na(drug)))
    finalCatalog = finalCatalog %>%
      bind_rows(curSubCatalog)
  }
  finalCatalog = finalCatalog %>%
    rename(Supplementary_Grading_Considerations = `Additional grading criteria`, Initial_Confidence_Grading = Initial, Final_Confidence_Grading = Final)
  manual_check_results = read_csv(paste0(NON_DATABASE_DIRECTORY, "/manual_check.csv"), guess_max = LARGE_NUMBER, show_col_types = FALSE) %>%
    select(drug, variant, Supplementary_Grading_Considerations, Final_Confidence_Grading) %>%
    mutate_at("Supplementary_Grading_Considerations", ~{str_replace(., "Additional grading evidence", "Evidence")})
  ## End of new block
  finalCatalog = applyManualChecks(finalCatalog, manual_check_results)
  write_csv(finalCatalog, file.path(OUTPUT_DIRECTORY, paste0("Final_graded_algorithm_catalogue_", format(Sys.Date(), "%d%b%Y"), "_withLoFs.csv")))
  reducedGradedCatalog = finalCatalog %>%
    select(c(drug, variant, Initial_Confidence_Grading, Supplementary_Grading_Considerations,	Final_Confidence_Grading))
  write_csv(reducedGradedCatalog, file = file.path(OUTPUT_DIRECTORY, paste0("List_of_graded_variants_", format(Sys.Date(), "%d%b%Y"), ".csv")))
  ## Additional processing: orphans and lineage stratification of all isolates with the mutation present
  if (augmentWithOrphansAndLineage) {
    orphanTab = getOrphanData(DATA_DIRECTORY, EXTRACTION_ID, minMAF = minMAF, minQ = minQ)
    lineageTab = getLineageData(DATA_DIRECTORY, EXTRACTION_ID, useSublineageData = TRUE)
    mainTab = fullDataset[["MAIN"]]
    fullTab = mainTab %>%
      bind_rows(orphanTab) %>%
      left_join(lineageTab, by = "sample_id") %>%
      pivot_wider(names_from = lineage, values_from = lineage, values_fill = list(lineage = 0), values_fn = length)
    colnames(fullTab)[-(1:ncol(mainTab))] = paste0("lineage_", colnames(fullTab)[-(1:ncol(mainTab))])
    fullCounts = fullTab %>%
      dplyr::filter(!het & !(sample_id %in% samplesToExclude$sample_id)) %>%
      group_by(drug, variant) %>%
      summarise(across(starts_with('lineage'), sum), 
                count = n(), orphan_count = sum(phenotype == "U"), .groups = "drop")
    finalCatalog = finalCatalog %>%
      full_join(fullCounts, by = c("drug", "variant")) %>%
      mutate(across(starts_with("lineage") | ends_with("count"), ~{ifelse(is.na(.), 0, .)}))
    cnames = colnames(finalCatalog)
    countNames = cnames[str_detect(cnames, "lineage_")]
    countNames = countNames[order(str_remove(countNames, "lineage_") %>% as.numeric())]
    cnames = c(setdiff(cnames, countNames), countNames)
    finalCatalog = finalCatalog %>%
      select(all_of(cnames))
    write_csv(finalCatalog, file.path(OUTPUT_DIRECTORY, paste0("Final_graded_algorithm_catalogue_", format(Sys.Date(), "%d%b%Y"), "_withLoFs_and_counts.csv")))
  }
  fullDataset
}

#' Extract orphan (no phenotype, genotype-only) data
#' @noRd
getOrphanData = function(DATA_DIRECTORY, EXTRACTION_ID, minMAF = MAF_THRESHOLD_REGULAR, minQ = QUALITY_THRESHOLD_STRICT) {
  orphanDir = file.path(str_remove(DATA_DIRECTORY, "/$"), str_remove(EXTRACTION_ID, "/$"), "orphan_genotypes")
  orphanData = read_csv(list.files(orphanDir, full.names = TRUE)[1], guess_max = LARGE_NUMBER, show_col_types = FALSE) %>%
    rename(drug = drug_name, gene = resolved_symbol, mutation = variant_category, effect = predicted_effect) %>%
    mutate(variant = paste0(gene, "_", mutation), phenotype = "U") %>%
    applyThresholds(minMAF = minMAF, lowMAFHet = TRUE, minQ = minQ, lowQHet = TRUE) %>%
    extractPositions(colName = "mutation", maxNumber = 2L) %>%
    computeRRDRInfo()
  orphanData
}

#' Extract per-isolate lineage data
#' @noRd
getLineageData = function(DATA_DIRECTORY, EXTRACTION_ID, useSublineageData = TRUE) {
  lineageCountsDir  = file.path(str_remove(DATA_DIRECTORY, "/$"), str_remove(EXTRACTION_ID, "/$"), "lineages_counts")
  lineageCountsFile = list.files(lineageCountsDir, full.names = TRUE)
  stopifnot(length(lineageCountsFile) == 1)
  lineageRaw  = read_csv(lineageCountsFile, col_types = cols(sample_id = col_double(), position = col_integer(), final_af = col_double(), .default = col_character()), show_col_types = FALSE)
  lineageData = computeDominantLineage(lineageRaw, mafThresholds = MAF_THRESHOLD_REGULAR, subLineage = useSublineageData)[[paste0("maf", MAF_THRESHOLD_REGULAR)]] %>%
    select(sample_id, dominant_lineage) %>%
    rename(lineage = dominant_lineage) %>%
    mutate_at("lineage", ~{ifelse(!(str_sub(., 1, 1) %in% as.character(1:6)), "Other", .)})
  lineageData
}

#' Run consistency checks on genotype data
#' @noRd
runGenotypeConsistencyTests = function(allGenotypes) {
  stopifnot(all(allGenotypes %>% dplyr::filter(gene ==  PZA_GENE) %>% pull(drug) == "Pyrazinamide"))
  stopifnot(all(allGenotypes %>% dplyr::filter(gene == RRDR_GENE) %>% pull(drug) == "Rifampicin"))
  ## the initial amino acid is always consistent for a given gene-position combination
  allGenotypes = allGenotypes %>%
    mutate(aa_start = ifelse(str_starts(mutation, "p.") & is.na(pos2), str_sub(mutation, 3, 5), NA))
  stopifnot(all(testConsistent(allGenotypes %>% dplyr::filter(!is.na(aa_start)), c("gene", "pos1"), consistentVars = c("aa_start"))[[1]]))
  ## the mutation is missing if and only if the effect is missing
  stopifnot(all((allGenotypes$mutation == "missing") == (allGenotypes$effect == "missing")))
  ## same effect and position for each variant-drug combination
  stopifnot(all(testConsistent(allGenotypes, c("drug", "variant"), consistentVars = c("effect", "position"))[[1]]))
  ## same tier for each gene-drug combination
  stopifnot(all(testConsistent(allGenotypes, c("drug", "gene"   ), consistentVars = "tier")[[1]]))
}

#' Run consistency checks on phenotype data
#' @noRd
runPhenotypeConsistencyTests = function(allPhenotypes) {
  ## Initial consistency checks on the phenotypes: phenotype consistent for sample-drug combinations per category
  stopifnot(testConsistent(allPhenotypes, c("sample_id", "category_phenotype", "drug"), consistentVars = "phenotype")[[1]])
  ## Initial consistency checks on the phenotypes: unique isolate-drug combinations between the WHO, ALL datasets
  stopifnot(anyDuplicated(allPhenotypes %>% dplyr::filter(category_phenotype %in% c("WHO", "ALL")) %>% select(drug, sample_id)) == 0)
}

#' Mark SOLO algorithm stage assignments in a result table
#' @noRd
markStages = function(inputTab) {
  inputTab = inputTab %>%
    mutate(stage = case_when(
      tier == 1 & !(effect %in% SILENT_EFFECTS) ~ 1L,
      tier == 2 & !(effect %in% SILENT_EFFECTS) ~ 2L,
      tier == 1 &  (effect %in% SILENT_EFFECTS) ~ 3L,
      .default                                  = 4L
    ))
  inputTab
}
