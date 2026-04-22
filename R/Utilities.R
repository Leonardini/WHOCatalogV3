#' Extract genotypes (geno = TRUE) or phenotypes (geno = FALSE) stored within inDir
#' @noRd
extractData = function(inDir, drugList = NULL, geno = FALSE) {
  ## Identify sub-directories
  useDirs = list.files(inDir)
  ## Only keep the sub-directories matching the drug list, if specified
  if (!is.null(drugList)) {
    useDirs %<>%
      intersect(paste0("drug_name=", drugList))
  }
  ## Update the drug list accordingly
  drugList = str_remove_all(useDirs, "drug_name=")
  ## Iteration counter
  N = length(useDirs)
  ## Initially empty output table
  outputTab = tibble()
  ## Iterate over drugs
  for (index in 1:N) {
    message(useDirs[[index]])
    drugDir = file.path(inDir, useDirs[[index]])
    if (geno) {
      ## For genotypes, list the tiers
      curLF = list.files(drugDir)
      ## Extract the tier identifiers, convert to integers
      curTiers = str_remove_all(curLF, "tier=") %>%
        as.integer()
      ## Inner iteration counter
      M = length(curTiers)
      for (ind in 1:M) {
        tierDir = file.path(drugDir, curLF[[ind]])
        ## Identify the file, check it is unique
        curFile = list.files(tierDir, full.names = TRUE)
        stopifnot(length(curFile) == 1)
        ## Parse it and append it to the output table, recording the drug and the tier
        outputTab %<>%
          bind_rows(read_csv(curFile, guess_max = LARGE_NUMBER, show_col_types = FALSE) %>%
                      mutate(drug = drugList[index], tier = curTiers[ind]))
      }
    } else {
      ## Identify the file, check it is unique
      curFile = list.files(drugDir, full.names = TRUE)
      stopifnot(length(curFile) == 1)
      ## Parse it and append it to the output table, recording the drug
      outputTab %<>%
        bind_rows(read_csv(curFile, guess_max = LARGE_NUMBER, show_col_types = FALSE) %>%
                    mutate(drug = drugList[index]))
    }
  }
  outputTab
}

#' Combine genotypes and phenotypes, then split them according to phenotype group
#' @noRd
mergeGenoPheno = function(Genotypes, Phenotypes, phenoGroups = NULL) {
  Phenotypes %<>%
    inner_join(PHENO_GROUPS, by = "category_phenotype")
  fullDataset = right_join(Genotypes, Phenotypes, by = c("sample_id", "drug"), relationship = "many-to-many") %>%
    mutate_at("het", ~{ ifelse(is.na(variant), TRUE, .) })
  fullDataset %<>%
    split(fullDataset$group)
  fullDataset
}

#' Test that all the groups defined by groupingVars agree on consistentVars in a Table
#' @noRd
testConsistent <- function(Table, groupingVars, consistentVars) {
  if (length(consistentVars) == 0) { return(TRUE) }
  stopifnot(all(groupingVars %in% colnames(Table)))
  stopifnot(all(consistentVars %in% colnames(Table)))
  checks   <- setNames(rep(FALSE, length(consistentVars)), consistentVars)
  problems <- setNames(vector("list", length(consistentVars)), consistentVars)
  # For each consistent variable, check for multiple unique values per group
  for (curVar in consistentVars) {
    problem_rows <- Table %>%
      group_by(across(all_of(groupingVars))) %>%
      summarise(n_unique = n_distinct(.data[[curVar]]), .groups = "drop") %>%
      dplyr::filter(n_unique > 1)
    if (nrow(problem_rows) == 0) {
      checks[curVar] <- TRUE
    } else {
      problems[[curVar]] <- problem_rows
    }
  }
  list(checks = checks, problems = problems)
}

#' Read an Excel file with two header rows into a data frame
#' If postfix = TRUE, appends the first header row (propagated rightwards) to the second header row.
#' If criticalString is not NULL, uses it to select the relevant parts of the first header row.
#' Otherwise, directly uses the second header row whenever it is non-missing.
#' @noRd
read_excel_2headers = function(inputFile, postfix = TRUE, criticalString = PREFIX_STRING) {
  # Read the full data, skipping the first one or two rows
  initTab = readxl::read_excel(inputFile, sheet = 1L, skip  = 3L, col_names = FALSE, .name_repair = "unique_quiet", guess_max = LARGE_NUMBER)
  # Read the first two rows (header rows) separately
  headers = readxl::read_excel(inputFile, sheet = 1L, n_max = 2L, col_names = FALSE, .name_repair = "unique_quiet")
  if (postfix) {
    # Propagate the first header row
    firstHeader <- headers %>% 
      slice(1) %>%
      unlist(use.names = FALSE)
    for (ind in 2:ncol(headers)) {
      if (is.na(firstHeader[ind])) {
        firstHeader[ind] = firstHeader[ind - 1]
      }
    }
    firstHeader[is.na(firstHeader)] = ""
    if (!is.null(criticalString)) {
      firstHeader[!str_starts(firstHeader, criticalString)] = ""
      firstHeader = str_remove_all(firstHeader, criticalString) %>%
        str_trim()
    }
  } else {
    firstHeader = rep("", ncol(headers))
  }
  # Second header row as character vector
  secondHeader = headers %>%
    slice(2) %>%
    unlist(use.names = FALSE)
  # Combine headers with "_"
  combinedHeader = paste(secondHeader, firstHeader, sep = "_") %>%
    str_remove_all("_$")
  # Assign the combined headers as column names
  colnames(initTab) = combinedHeader
  initTab
}

#' Recode values in initTab according to the additionally supplied two-column manualTab
#' @noRd
recodeValues = function(initTab, manualTab) {
  stopifnot(ncol(manualTab) == 2)
  coln = colnames(manualTab)
  stopifnot(all(str_ends(coln, "_first", negate = TRUE)) && all(str_ends(coln, "_second", negate = TRUE)))
  ## Create the completed version of the sub-table, taking the merged column's values from manualTab
  otherTab = initTab %>%
    inner_join(manualTab, by = coln[1], suffix = c("_first", "_second")) %>%
    rename_with(.fn = function(x) { str_remove_all(x, "\\_second$") }) %>%
    select(!any_of(paste0(coln[2], "_first")))
  ## Now remove the problematic entries in the original table and add them back in from the completed sub-table
  initTab %<>%
    anti_join(otherTab, by = coln[1]) %>%
    bind_rows(otherTab)
  initTab
}

#' Simplify duplicate columns produced by a merge
#' The second column is used when the first is NA. If warn = TRUE, a warning is printed when any
#' pair of parallel columns has non-missing elements in common. If add = TRUE, numeric column pairs
#' are summed; non-numeric pairs are concatenated with a semicolon separator.
#' @noRd
adjustDuplicateColumns = function(initTab, suffixes = c(".x", ".y"), warn = TRUE, add = FALSE) {
  coln = colnames(initTab)
  coln1 = coln[str_ends(coln, suffixes[1])]
  coln2 = coln[str_ends(coln, suffixes[2])]
  coln1Short = str_sub(coln1, end = -(nchar(suffixes[1]) + 1))
  coln2Short = str_sub(coln2, end = -(nchar(suffixes[2]) + 1))
  stopifnot(all(sort(coln1Short) == sort(coln2Short)))
  for (ind in seq_along(coln1Short)) {
    Col = coln1Short[ind]
    vec1 = initTab %>%
      pull(all_of(paste0(Col, suffixes[1])))
    vec2 = initTab %>%
      pull(all_of(paste0(Col, suffixes[2])))
    if (warn && any(vec1 != vec2, na.rm = TRUE)) {
      warning(paste("column", Col, "has", sum(vec1 != vec2, na.rm = TRUE), "conflicting entries!"))
    }
    if (add) {
      if (is.numeric(vec1) && is.numeric(vec2)) {
        vec12 = ifelse(is.na(vec1), 0, vec1) + ifelse(is.na(vec2), 0, vec2)
        ### vec12[is.na(vec1) & is.na(vec2)] = NA
      } else {
        vec12 = vec1
        vec12[!is.na(vec2)] = vec2[!is.na(vec2)]
        vec12[!is.na(vec1) & !is.na(vec2)] = paste(vec1[!is.na(vec1) & !is.na(vec2)], vec2[!is.na(vec1) & !is.na(vec2)], sep = "; ")
      }
    } else {
      vec12 = ifelse(!is.na(vec1), vec1, vec2)
    }
    initTab %<>%
      bind_cols(enframe(vec12, name = NULL, value = Col))
  }
  initTab %<>%
    select(!any_of(c(coln1, coln2)))
  initTab
}

#' @noRd
mergeCountColumns = function(initTab) adjustDuplicateColumns(initTab, warn = FALSE, add = TRUE)

#' Prepare a masked dataset according to four criteria: synonymous, tier 2, neutral, and Pool
#' Silent = TRUE masks silent variants; Tier2 = TRUE masks tier 2 variants; Neutral = TRUE masks
#' neutral variants; Pool can be NA or a name masking corresponding poolable variants. SOnly = TRUE
#' marks rather than removes variants occurring only in S samples, and removes samples with hets.
#' @noRd
prepMask = function(inputTab, Silent = FALSE, Tier2 = TRUE, Neutral = TRUE, Pool = NA, SOnly = TRUE) {
  if (Silent) {
    inputTab %<>%
      dplyr::filter(effect == "missing" | !(effect %in% SILENT_EFFECTS))
  }
  if (Tier2) {
    inputTab %<>%
      dplyr::filter(tier != 2)
  }
  if (Neutral) {
    inputTab %<>%
      dplyr::filter(!neutral)
  }
  if (!is.na(Pool)) { ## Minor fix for a bug identified by Sacha, 3rd July 2024
    inputTab %<>%
      dplyr::filter(str_ends(variant, Pool) | effect == "missing" | (effect %in% POOLED_EFFECTS[[Pool]] & het) | !(effect %in% POOLED_EFFECTS[[Pool]]))
  }
  if (SOnly) {
    inputTab %<>%
      group_by(drug, variant) %>%
      mutate(SOnly = (sum(!het & phenotype == "S") == sum(!het))) %>%
      ungroup() %>%
      removeSamplesWithHets() %>%
      select(drug, sample_id, variant, phenotype, SOnly) %>%
      distinct()
  } else {
    inputTab %<>%
      select(drug, sample_id, variant, phenotype, het) %>%
      distinct()
  }
  inputTab
}

#' Identify sample-drug pairs that are not explained by the current variant classification
#' @noRd
findUnexplained = function(inputTab) {
  inputTab %<>%
    group_by(sample_id, drug) %>%
    mutate(satisfying = all(phenotype == "S") || (all(is.na(class) | class == "S"))) %>%
    ungroup() %>%
    dplyr::filter(satisfying) %>%
    select(-satisfying)
  inputTab
}

#' Remove sample-drug pairs that contain at least one heterozygous call
#' @noRd
removeSamplesWithHets = function(inputTab) {
  inputTab %<>%
    group_by(sample_id, drug) %>%
    mutate(satisfying = (all(!het))) %>%
    ungroup() %>%
    dplyr::filter(satisfying) %>%
    select(-satisfying)
  inputTab
}

#' Convert a vector to its logical form, treating NA as FALSE
#' @noRd
convertToLogical = function(x) {
  output = ifelse(is.na(x), FALSE, as.logical(x))
  output
}

#' Expand list columns of binom.test results into point estimate and CI bound columns
#' @noRd
expandBinomCIs = function(df) {
  df %>%
    mutate(across(where(is.list), list(lb = ~{map_dbl(., ~{.$conf.int[1]})}, ub = ~{map_dbl(., ~{.$conf.int[2]})}))) %>%
    mutate(across(where(is.list), ~{map_dbl(., ~{.$estimate})}))
}

#' Compute binomial confidence intervals for each distinct pair of numerator and denominator counts
#' @noRd
binomLookup = function(df, num_col, denom_col, stat_name) {
  df %>%
    distinct(across(all_of(c(num_col, denom_col))), .keep_all = FALSE) %>%
    mutate(!!stat_name := map2(.data[[num_col]], .data[[num_col]] + .data[[denom_col]], safeBinomTest))
}

#' Extract up to maxNumber integer positions from a specified column in a data frame
#' @param inputTable A data frame containing the column to extract positions from.
#' @param colName Name of the column containing position strings.
#' @param maxNumber Maximum number of integer positions to extract per row.
#' @noRd
extractPositions = function(inputTable, colName = "mutation", maxNumber = 2L) {
  stopifnot(is.data.frame(inputTable))
  stopifnot(is.character(colName) && length(colName) == 1)
  stopifnot(colName %in% names(inputTable))
  stopifnot(is.numeric(maxNumber) && maxNumber > 0)
  # Extract all integers from the target column
  extracted <- str_extract_all(inputTable[[colName]], "-?\\d+")
  # Check for rows with too many matches
  num_matches <- map_int(extracted, length)
  if (any(num_matches > maxNumber)) {
    offending_rows <- which(num_matches > maxNumber)
    stop(paste0("Rows ", paste(offending_rows, collapse = ", "), " contain more than ", maxNumber, " integers."))
  }
  # Create pos1 ... posK columns, where K = maxNumber
  for (i in seq_len(maxNumber)) {
    inputTable[[paste0("pos", i)]] <- map_int(extracted, ~ {
      if (length(.x) >= i) as.integer(.x[i]) else NA_integer_
    })
  }
  inputTable
}
