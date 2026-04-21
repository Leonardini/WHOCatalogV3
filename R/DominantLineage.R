#' Extract the sub-lineage from a lineage numbering string (up to two dot-separated levels, "lineage" prefix removed)
#' @noRd
toSubLineage = function(s) {
  ifelse(is.na(s) | s == "", NA_character_,
    str_remove(str_extract(s, "^[^.]+(?:\\.[^.]+)?"), "^lineage"))
}

#' Summarise marker data for a single MAF threshold, returning dominant sub-lineage per sample
#' @noRd
summarizeByMAF = function(inputTab, maf, subLineage) {
  allSamples = tibble(sample_id = unique(inputTab$sample_id))
  kept = inputTab %>%
    dplyr::filter(final_af >= maf) %>%
    select(sample_id, lineage_numbering) %>%
    mutate(lineage_level = if (subLineage) toSubLineage(lineage_numbering) else
      str_remove(str_extract(lineage_numbering, "^[^.]+"), "^lineage"))
  if (nrow(kept) == 0) {
    return(allSamples %>%
      mutate(dominant_lineage = "", ties = "") %>%
      arrange(sample_id))
  }
  counts = kept %>%
    group_by(sample_id, lineage_level) %>%
    summarise(n = n(), .groups = "drop")
  top = counts %>%
    group_by(sample_id) %>%
    dplyr::filter(n == max(n)) %>%
    summarise(all_lineages = list(sort(lineage_level)), .groups = "drop") %>%
    mutate(
      dominant_lineage = map_chr(all_lineages, ~ .x[[1]]),
      ties             = map_chr(all_lineages, ~ paste(.x[-1], collapse = ","))
    ) %>%
    select(-all_lineages)
  allSamples %>%
    left_join(top, by = "sample_id") %>%
    mutate(
      dominant_lineage = replace_na(dominant_lineage, ""),
      ties             = replace_na(ties, "")
    ) %>%
    arrange(sample_id)
}

#' Compute the dominant sub-lineage for each sample from lineage marker data
#'
#' For each MAF threshold, rows with \code{final_af >= threshold} are retained, the
#' \code{lineage_numbering} field is collapsed to the sub-lineage level (two
#' dot-separated segments with any leading \code{"lineage"} prefix stripped), and
#' the lineage supported by the most markers is chosen as the dominant one.  When
#' multiple lineages tie, the alphabetically first is reported as
#' \code{dominant_lineage} and the rest appear in \code{ties} (comma-separated).
#' Samples with no rows passing the threshold are still included with empty strings.
#'
#' @param inputTab A data frame with columns \code{sample_id}, \code{lineage_numbering}, and \code{final_af}.
#' @param mafThresholds A numeric vector of MAF thresholds (default \code{c(0.25, 0.75, 0.9)}).
#' @param subLineage Logical; if \code{TRUE} (default) extracts the sub-lineage (up to two dot-separated levels); if \code{FALSE} extracts only the broad lineage (one level).
#' @param outDir Optional directory for writing one output CSV per threshold; \code{NULL} skips file writing.
#' @return A named list of data frames (names like \code{"maf0.25"}), each with columns \code{sample_id}, \code{dominant_lineage}, and \code{ties}.
#' @export
computeDominantLineage = function(inputTab, mafThresholds = c(0.25, 0.75, 0.9),
                                  subLineage = TRUE, outDir = NULL) {
  stopifnot(is.data.frame(inputTab))
  missing_cols = setdiff(c("sample_id", "lineage_numbering", "final_af"), colnames(inputTab))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(sort(missing_cols), collapse = ", ")))
  }
  results = setNames(
    lapply(mafThresholds, function(maf) summarizeByMAF(inputTab, maf, subLineage)),
    paste0("maf", mafThresholds)
  )
  if (!is.null(outDir)) {
    for (maf in mafThresholds) {
      write_csv(results[[paste0("maf", maf)]],
                file.path(outDir, paste0("dominant_lineage_maf", maf, ".csv")))
    }
  }
  results
}
