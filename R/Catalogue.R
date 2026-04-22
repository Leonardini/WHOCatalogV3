#' Apply a catalogue file to an input data table, replacing grades and applying RRDR/LoF rules
#' @param inputData A data frame with columns drug, effect, gene, mutation, and variant.
#' @param catalogueFile Path to a catalogue file (.csv or .xlsx).
#' @param version Catalogue version to apply; controls DLM/PMD cross-resistance handling.
#' @inheritParams mainDriver
#' @return \code{inputData} with additional columns: \code{het} (logical, TRUE if the variant fails
#'   MAF or quality thresholds), \code{RRDR_NON_SILENT} (logical), \code{LoF_candidate} (logical,
#'   TRUE for novel LoF mutations in genes that have a graded pooled LoF), and \code{Final} (integer
#'   catalogue grade, NA if the variant is absent from the catalogue).
#' @export
applyCatalogue = function(inputData, catalogueFile, minMAF = MAF_THRESHOLD_REGULAR, lowMAFHet = TRUE, 
                          minQ = QUALITY_THRESHOLD_STRICT, lowQHet = TRUE,
                          LoF = TRUE, version = CURR_VERSION) {
  # Check if the input data has the required columns
  if (!all(c("drug", "effect", "gene", "mutation", "variant") %in% colnames(inputData))) {
    stop("Input data must contain 'drug', 'effect', 'gene', 'mutation' and 'variant' columns.")
  }
  # Compute the het variable on the input data using the MAF and quality cutoffs
  inputData = inputData %>% 
    applyThresholds(minMAF = minMAF, lowMAFHet = lowMAFHet, minQ = minQ, lowQHet = lowQHet) %>%
    mutate(across(starts_with("het"), ~{replace_na(., TRUE)}))
  # Compute the RRDR_NON_SILENT variable if needed
  inputData = inputData %>%
    computeRRDRInfo()
  # Load the catalogue
  if (str_ends(catalogueFile, 'xlsx')) {
    catalogueData = read_excel_2headers(catalogueFile, postfix = TRUE)
    coln = unique(colnames(catalogueData))
    catalogueData = catalogueData[, coln] %>%
      as_tibble() %>%
      mutate(Final = match(`FINAL CONFIDENCE GRADING`, GRADES))
  } else {
    catalogueData = read_csv(catalogueFile, guess_max = LARGE_NUMBER, show_col_types = FALSE) %>%
      mutate(Final = match(Final_Confidence_Grading, GRADES))
  }
  catalogueData = catalogueData %>%
    select(drug, variant, Final)
  ### For v2, DLM/PMD cross-resistance, add SPECIAL_DLM_VAR to PMD graded 2
  if (version == PREV_VERSION) {
    catalogueData = catalogueData %>%
      bind_rows(tibble(drug = "Pretomanid", variant = SPECIAL_DLM_VAR, Final = 2L))
  }
  if (LoF) {
    LoFTab = catalogueData %>%
      dplyr::filter(str_detect(variant, "_LoF")) %>%
      mutate(gene = str_remove(variant, "_LoF")) %>%
      dplyr::filter(Final <= RESISTANCE_GRADE_MAX)
    ### For v2, DLM/PMD cross-resistance, add genes that have an LoF graded 1/2
    if (version == PREV_VERSION) {
      extraTab = LoFTab %>%
        dplyr::filter(drug == "Delamanid") %>%
        mutate(drug = "Pretomanid")
      LoFTab = bind_rows(LoFTab, extraTab)
    }
  }
  # Merge the input data with the catalogue and apply the expert rule on RRDR
  mergedData = left_join(inputData, catalogueData, by = c("drug", "variant")) %>%
    mutate(Final = ifelse(RRDR_NON_SILENT & is.na(Final), 2L, Final))
  # Apply the expert rule on LOFs: if a pooled LoF is grade <= 2, any NOVEL LoF mutation in the same drug-gene combination is grade 2
  if (LoF) {
    mergedData = mergedData %>%
      mutate(LoF_candidate = (paste0(gene, drug, sep = "_") %in% paste0(LoFTab$gene, LoFTab$drug, sep = "_") 
                              & effect %in% POOLED_EFFECTS[[LOF_LABEL]] & is.na(Final))) %>%
      mutate(Final = ifelse(LoF_candidate, 2L, Final))
  }
  # Return the merged data
  return(mergedData)
}

#' Apply MAF and quality thresholds to input data
#' If lowMAFHet = TRUE, low-MAF variants are treated as hets; otherwise they are filtered out.
#' If lowQHet = TRUE, low-quality variants are treated as hets; otherwise they are filtered out.
#' @noRd
applyThresholds = function(inputData, minMAF = MAF_THRESHOLD_REGULAR, lowMAFHet = TRUE, 
                           minQ = QUALITY_THRESHOLD_STRICT, lowQHet = TRUE) {
  if ('max(af)' %in% colnames(inputData) && !is.na(minMAF)) {
    if (lowMAFHet) {
      inputData = inputData %>%
        mutate(het = (`max(af)` < minMAF | variant == MISSING_VARIANT))
    } else {
      inputData = inputData %>%
        dplyr::filter(is.na(`max(af)`) | `max(af)` >= minMAF)
    }
  } else {
    warning('no max(af) column found in input data; skipping MAF filtering.')
  }
  if ('max(quality)' %in% colnames(inputData) && !is.na(minQ)) {
    if (lowQHet) {
      inputData = inputData %>%
        mutate(het = het | (`max(quality)` < minQ))
    } else {
      inputData = inputData %>%
        dplyr::filter(is.na(`max(quality)`) | `max(quality)` >= minQ)
    }
  } else {
    warning('no max(quality) column found in input data; skipping quality filtering.')
  }
  return(inputData)
}

#' Apply MAF/quality thresholds, extract positions, and compute RRDR non-silent status
#' @noRd
applyVariantPostprocessing = function(inputData, minMAF = MAF_THRESHOLD_REGULAR, lowMAFHet = TRUE, minQ = QUALITY_THRESHOLD_STRICT, lowQHet = TRUE) {
  inputData %>%
    applyThresholds(minMAF = minMAF, lowMAFHet = lowMAFHet, minQ = minQ, lowQHet = lowQHet) %>%
    extractPositions(colName = "mutation", maxNumber = 2L) %>%
    computeRRDRInfo()
}

# This function computes the RRDR_NON_SILENT variable if necessary
#' Compute RRDR (rifampicin resistance-determining region) information for each variant
#' @noRd
computeRRDRInfo = function(inputData) {
  if (!"RRDR_NON_SILENT" %in% colnames(inputData)) {
    if (!("pos1" %in% colnames(inputData) && "pos2" %in% colnames(inputData))) {
      inputData = inputData %>%
        extractPositions(colName = "mutation", maxNumber = 2L)
    }
    inputData = inputData %>%
      mutate(RRDR_NON_SILENT = (gene == RRDR_GENE & (pos1 %in% RRDR_INT | pos2 %in% RRDR_INT) & !(effect %in% SILENT_EFFECTS)))
  }
  return(inputData)
}
