## source("SOLOport/SensSpec.R")
# Load required libraries
library(arrow)
library(httr)
library(xml2)
library(jsonlite)
library(dplyr)
library(purrr)
library(magrittr)
library(tidyverse)
library(readxl)

conflicts_prefer(dplyr::rename)
conflicts_prefer(dplyr::slice)

ACCESSIONS = c(
  PROJECT      = "^PRJ(E|D|N)[A-Z][0-9]+",
  STUDY        = "^(E|D|S)RP[0-9]{6,}",
  BIOSAMPLE    = "^SAM(D|N|EA)[0-9]+",
  SAMPLE       = "^(E|D|S)RS[0-9]{6,}",
  EXPERIMENT   = "^(E|D|S)RX[0-9]{6,}",
  RUN          = "^(E|D|S)RR[0-9]{6,}",
  ANALYSIS     = "^(E|D|S)RZ[0-9]{6,}",
  ASSEMBLY     = "^GC(A|F)[\\_][0-9]{9}[\\.]?[0-9]*",
  PROTEIN      = "^[A-Z]{3}[0-9]{5,7}[\\.][0-9]+",
  SCAFFOLD     = "^[A-Z]{1,2}[0-9]{5,6}[\\.]?[0-9]*",
  CONTIG       = "^[A-Z]{4,6}[0-9]{2}[S]?[0-9]{6,9}"
)

EXTRACTION_ID = "2025-12-16T22-09-33.644305_jr_639ca1c8b567dc6a6e49b9097bd7c369c74bb0ba0a5749c1d9e886890abaa333"
DATA_DIRECTORY = "SOLO Algorithm Input files/DATABASE EXTRACTION files/"

EXTRA_ID = "2025-09-21T17-50-00.187620_jr_a7b222b6f60a7318d8214906feeda9bbda0f9a97ac365f46af63242889b3dbb1"

analyseIDs = function(checkPrevious = FALSE, minLength = 5L, keepInfo = FALSE) {
  setwd(paste0(DATA_DIRECTORY, EXTRA_ID))
  setwd("identifiers")
  Tab = read_csv(list.files(pattern = "run")[1], show_col_types = FALSE, guess_max = 1e5)
  if (checkPrevious) { ## Check that we have all the IDs in the previous dataset
    prevData = read_csv(paste0("../../../../Results/", EXTRACTION_ID, "/completeDataset.csv"),
                        show_col_types = FALSE, guess_max = 1e5)
    UIDs = sort(unique(prevData$sample_id))
    stopifnot(all(UIDs %in% Tab$sample_id))
  }
  ## Change column names and do some parsing
  colnames(Tab)[-1] = c("name", "fastq", "library")
  ## Split regular names into sample names, run names, and other names
  allNames = Tab$name
  allNames = str_split(allNames, ", ")
  allNames = lapply(allNames, function(x) {str_split(x, "__") %>% unlist() %>% unique() %>% sort() })
  sampleNames = lapply(allNames, function(x) {
    x[str_detect(x, ACCESSIONS[["SAMPLE"]]) | str_detect(x, ACCESSIONS[["BIOSAMPLE"]])] })
  runNames = lapply(allNames, function(x) { x[str_detect(x, ACCESSIONS[["RUN"]])] })
  otherNames  = lapply(1:length(allNames), function(x) { setdiff(allNames[[x]], c(sampleNames[[x]], runNames[[x]])) })
  ## Split fastq names into R1 and R2 and collect unique values
  fastqNames = Tab$fastq %>%
    lapply(function(x) { str_split(x, "; ") %>% unlist() %>% str_remove_all("_R[0-9].fastq.gz") %>% unique() %>% sort() })
  stopifnot(max(lengths(fastqNames)) <= 1)
  frunNames = lapply(fastqNames, function(x) { x[str_detect(x, ACCESSIONS[["RUN"]])] })
  fotherNames = lapply(1:length(fastqNames), function(x) { setdiff(fastqNames[[x]], frunNames[[x]]) })
  ## Split library names into run names and other names
  libNames = Tab$library %>%
    lapply(function(x) { str_split(x, ", ") %>% unlist() %>% unique() %>% sort() })
  xrunNames = lapply(libNames, function(x) { x[str_detect(x, ACCESSIONS[["RUN"]])] })
  xotherNames = lapply(1:length(libNames), function(x) { setdiff(libNames[[x]], xrunNames[[x]]) })
  ## Collect run names and other names together
  allSampleNames = sampleNames
  allRunNames = mapply(function(x, y, z) { sort(unique(c(x, y, z))) }, runNames, frunNames, xrunNames, SIMPLIFY = FALSE)
  allOtherNames = mapply(function(x, y, z) { sort(unique(c(x, y, z))) }, otherNames, fotherNames, xotherNames, SIMPLIFY = FALSE)
  Tab$name = allOtherNames
  Tab$fastq = allRunNames
  Tab$library = allSampleNames
  fullTab = Tab %>%
    unnest(name, keep_empty = TRUE) %>%
    unnest(fastq, keep_empty = TRUE) %>%
    unnest(library, keep_empty = TRUE)
  ## Check that we have all the identifiers from Maha's file:
  allIDs = tibble(sample_id = rep(fullTab$sample_id, 3), name = trimws(c(fullTab$library, fullTab$fastq, fullTab$name))) %>%
    filter(!is.na(name) & name != "") %>%
    distinct()
  problemIDs = allIDs %>%
    group_by(name) %>% 
    mutate(N  = n()) %>% 
    ungroup() %>% 
    filter(N > 1) %>%
    arrange(nchar(name), name)
  write_csv(problemIDs, file = "problematicIDs.csv")
  mahaFile = read_xlsx("../../../../Maha Regression Supplementary File.xlsx") %>%
    select(Sample_ID, BioSample, Sample_Name, Lineage) %>%
    mutate_all(trimws) %>%
    mutate_at("Sample_ID", as.numeric)
  mahaIDs = tibble(sample_id = rep(mahaFile$Sample_ID, 2), name = trimws(c(mahaFile$BioSample, mahaFile$Sample_Name))) %>%
    filter(!is.na(name) & name != "") %>%
    mutate(short_name = paste0(sapply(str_split(name, "_"), dplyr::nth, -3), "_", sapply(str_split(name, "_"), dplyr::nth, -2))) %>%
    mutate(short_name = toupper(ifelse(str_count(name, "_") >= 2, short_name, name))) %>%
    mutate(minimal_name = toupper(ifelse(str_detect(short_name, "_"), sapply(str_split(short_name, "_"), dplyr::last), name)))
  mahaProblems = mahaIDs %>%
    distinct() %>%
    group_by(name) %>% 
    mutate(N  = n()) %>% 
    ungroup() %>% 
    filter(N > 1) %>%
    arrange(nchar(name), name)
  stopifnot(nrow(mahaProblems) == 0)
  mahaFile$status = "unmatched"
  initialMatch = allIDs %>%
    inner_join(mahaIDs, by = "name", relationship = "many-to-many") %>%
    dplyr::rename(v3_sample_id = sample_id.x, v2_sample_id = sample_id.y) %>%
    distinct(v2_sample_id, v3_sample_id, .keep_all = TRUE) %>%
    select(v2_sample_id, v3_sample_id, name)
  mahaFile$status[mahaFile$Sample_ID %in% initialMatch$v2_sample_id] = "matched by full name"
  minimalMatch = allIDs %>%
    filter(nchar(name) >= minLength) %>%
    inner_join(mahaIDs %>% filter(minimal_name != name), by = join_by(name == minimal_name), relationship = "many-to-many") %>%
    dplyr::rename(v3_sample_id = sample_id.x, v2_sample_id = sample_id.y) %>%
    distinct(v2_sample_id, v3_sample_id, .keep_all = TRUE) %>%
    select(v2_sample_id, v3_sample_id, name)
  initialMatch %<>%
    rbind(minimalMatch)
  mahaFile$status[mahaFile$Sample_ID %in% minimalMatch$v2_sample_id] = "matched by short name"
  SRA_DIR = "~/Downloads/FINDConsulting/SRATables/"
  initDir = getwd()
  setwd(SRA_DIR)
  mahaUnmatched = mahaIDs %>%
    filter(!(sample_id %in% initialMatch$v2_sample_id)) %>%
    arrange(sample_id, nchar(name), name) %>%
    distinct()
  mahaUnmatchedBio = mahaUnmatched %>%
    filter(str_detect(name, ACCESSIONS[["BIOSAMPLE"]]))
  mahaUnmatched = mahaUnmatched %>%
    filter(!(sample_id %in% mahaUnmatchedBio$sample_id)) %>%
    arrange(nchar(name), name)

  summaryTabV3 = allIDs %>%
    filter(!(sample_id %in% initialMatch$v3_sample_id)) %>%
    mutate(prefix = str_extract(name, "[A-Za-z]+")) %>%
    group_by(prefix) %>%
    mutate(N = n_distinct(sample_id)) %>%
    filter(N > 30) %>%
    slice(1) %>%
    ungroup() %>%
    rename(example = name) %>%
    filter(!is.na(prefix)) %>%
    select(prefix, N, example, sample_id) %>%
    arrange(-N)
  
  summaryTabV2 = mahaUnmatched %>%
    filter(sample_id %in% (mahaFile %>% filter(status == "unmatched") %>% pull(Sample_ID))) %>%
    mutate(prefix = str_extract(short_name, "[A-Za-z]+")) %>%
    group_by(prefix) %>%
    mutate(N = n()) %>%
    filter(N > 30) %>%
    slice(1) %>%
    ungroup() %>%
    rename(example = short_name) %>%
    filter(!is.na(prefix)) %>%
    select(prefix, N, example, sample_id) %>%
    arrange(-N)
  if (keepInfo) {
    allInfo = tibble()
  }
  Tab1 = read_csv("SraRunTableCPATH.csv", show_col_types = FALSE, guess_max = 1e5) %>%
    mutate_all(trimws) %>%
    select(Run, BioSample, Experiment, `Sample Name`) %>%
    mutate(Sample_Name = paste0("SPE_", `Sample Name`)) %>%
    select(-`Sample Name`)
  extra1 = mahaUnmatched %>%
    inner_join(Tab1, by = join_by(name == Sample_Name), relationship = "many-to-many")
  mahaFile$status[mahaFile$Sample_ID %in% extra1$sample_id] = "found in CPATH"
  if (keepInfo) {
    allInfo %<>%
      rbind(extra1)
  }
  extra1 = extra1 %>%
    select(sample_id, Run) %>%
    rename(v2_sample_id = sample_id) %>%
    inner_join(allIDs, by = join_by(Run == name), relationship = "many-to-many") %>%
    rename(v3_sample_id = sample_id, name = Run) %>%
    distinct(v2_sample_id, v3_sample_id, .keep_all = TRUE) %>%
    select(v2_sample_id, v3_sample_id, name)
  initialMatch %<>%
    rbind(extra1)
  mahaFile$status[mahaFile$Sample_ID %in% extra1$v2_sample_id] = "matched after CPATH"
  
  Tab2 = read_csv("SraRunTableCRYPTIC.csv", show_col_types = FALSE, guess_max = 1e5) %>%
    mutate_all(trimws) %>%
    select(Run, BioSample, Experiment, `Library Name`) %>%
    mutate(Sample_Name = sapply(str_split(`Library Name`, " "), dplyr::last) %>% toupper()) %>%
    select(-`Library Name`) %>%
    filter(str_ends(Sample_Name, '-1-1')) %>%
    mutate_at("Sample_Name", ~{str_remove(., '-1-1$')})
  extra2 = mahaUnmatched %>%
    inner_join(Tab2, by = join_by(short_name == Sample_Name)) 
  mahaFile$status[mahaFile$Sample_ID %in% extra2$sample_id] = "found in CRYPTIC"
  if (keepInfo) {
    allInfo %<>%
      rbind(extra2)
  }
  extra2 = extra2 %>%
    select(sample_id, Run) %>%
    rename(v2_sample_id = sample_id) %>%
    inner_join(allIDs, by = join_by(Run == name), relationship = "many-to-many") %>%
    rename(v3_sample_id = sample_id, name = Run) %>%
    distinct(v2_sample_id, v3_sample_id, .keep_all = TRUE) %>%
    select(v2_sample_id, v3_sample_id, name)
  initialMatch %<>%
    rbind(extra2)
  mahaFile$status[mahaFile$Sample_ID %in% extra2$v2_sample_id] = "matched after CRYPTIC"
  
  Tab3 = read_csv("SraRunTableBORSTEL.csv", show_col_types = FALSE, guess_max = 1e5) %>%
    mutate_all(trimws) %>%
    select(Run, BioSample, Experiment, `Library Name`, Sample_name) %>%
    mutate(Sample_name = toupper(paste0(Sample_name, "_", `Library Name`))) %>%
    select(-`Library Name`)
  extra3 = mahaUnmatched %>%
    inner_join(Tab3, by = join_by(short_name == Sample_name)) 
  mahaFile$status[mahaFile$Sample_ID %in% extra3$sample_id] = "found in BORSTEL"
  if (keepInfo) {
    allInfo %<>%
      rbind(extra3)
  }
  extra3 = extra3 %>%
    select(sample_id, Run) %>%
    rename(v2_sample_id = sample_id) %>%
    inner_join(allIDs, by = join_by(Run == name), relationship = "many-to-many") %>%
    rename(v3_sample_id = sample_id, name = Run) %>%
    distinct(v2_sample_id, v3_sample_id, .keep_all = TRUE) %>%
    select(v2_sample_id, v3_sample_id, name)
  initialMatch %<>%
    rbind(extra3)
  mahaFile$status[mahaFile$Sample_ID %in% extra3$v2_sample_id] = "matched after BORSTEL"
  
  Tab4 = read_csv("SraRunTablePHO.csv", show_col_types = FALSE, guess_max = 1e5) %>%
    mutate_all(trimws) %>%
    select(Run, BioSample, Experiment, `Library Name`) %>%
    mutate(Sample_Name = paste0(sapply(str_split(`Library Name`, ". "), dplyr::first), "_", sapply(str_split(`Library Name`, ". "), dplyr::nth, -2))) %>%
    select(-`Library Name`) %>%
    mutate_at("Sample_Name", toupper)
  extra4 = mahaUnmatched %>%
    inner_join(Tab4, by = join_by(short_name == Sample_Name)) 
  mahaFile$status[mahaFile$Sample_ID %in% extra4$sample_id] = "found in PHO"
  if (keepInfo) {
    allInfo %<>%
      rbind(extra4)
  }
  extra4 = extra4 %>%
    select(sample_id, Run) %>%
    rename(v2_sample_id = sample_id) %>%
    inner_join(allIDs, by = join_by(Run == name), relationship = "many-to-many") %>%
    rename(v3_sample_id = sample_id, name = Run) %>%
    distinct(v2_sample_id, v3_sample_id, .keep_all = TRUE) %>%
    select(v2_sample_id, v3_sample_id, name)
  initialMatch %<>%
    rbind(extra4)
  mahaFile$status[mahaFile$Sample_ID %in% extra4$v2_sample_id] = "matched after PHO"
  
  Tab5p = read_csv("SraRunTableMONASH+.csv", show_col_types = FALSE, guess_max = 1e5) %>%
    select(Run, BioSample, Experiment, `Library Name`)
  Tab5 = read_csv("SraRunTableMONASH.csv", show_col_types = FALSE, guess_max = 1e5) %>%
    select(Run, BioSample, Experiment, `Library Name`) %>%
    bind_rows(Tab5p) %>%
    mutate_all(trimws) %>%
    mutate(Sample_Name = paste0(sapply(str_split(`Library Name`, ". "), dplyr::first), "_", sapply(str_split(`Library Name`, ". "), dplyr::nth, -2))) %>%
    select(-`Library Name`) %>%
    mutate_at("Sample_Name", toupper)
  extra5 = mahaUnmatched %>%
    inner_join(Tab5, by = join_by(short_name == Sample_Name)) 
  mahaFile$status[mahaFile$Sample_ID %in% extra5$sample_id] = "found in MONASH"
  if (keepInfo) {
    allInfo %<>%
      rbind(extra5)
  }
  extra5 = extra5 %>%
    select(sample_id, Run) %>%
    rename(v2_sample_id = sample_id) %>%
    inner_join(allIDs, by = join_by(Run == name), relationship = "many-to-many") %>%
    rename(v3_sample_id = sample_id, name = Run) %>%
    distinct(v2_sample_id, v3_sample_id, .keep_all = TRUE) %>%
    select(v2_sample_id, v3_sample_id, name)
  initialMatch %<>%
    rbind(extra5)
  mahaFile$status[mahaFile$Sample_ID %in% extra5$v2_sample_id] = "matched after MONASH"
  
  Tab6m = read_csv("SraRunTableOSR-.csv", show_col_types = FALSE, guess_max = 1e5) %>%
    select(Run, BioSample, Experiment, `Library Name`)
  Tab6p = read_csv("SraRunTableOSR+.csv", show_col_types = FALSE, guess_max = 1e5) %>%
    select(Run, BioSample, Experiment, `Library Name`)
  Tab6 = read_csv("SraRunTableOSR.csv", show_col_types = FALSE, guess_max = 1e5) %>%
    select(Run, BioSample, Experiment, `Library Name`) %>%
    bind_rows(Tab6p) %>%
    bind_rows(Tab6m) %>%
    mutate_all(trimws) %>%
    mutate(Sample_Name = paste0(sapply(str_split(`Library Name`, ". "), dplyr::nth, -2), "_", sapply(str_split(`Library Name`, ". "), dplyr::first))) %>%
    select(-`Library Name`) %>%
    mutate_at("Sample_Name", ~{str_remove_all(., "_06TB") %>% toupper})
  extra6 = mahaUnmatched %>%
    inner_join(Tab6, by = join_by(short_name == Sample_Name)) 
  mahaFile$status[mahaFile$Sample_ID %in% extra6$sample_id] = "found in OSR"
  if (keepInfo) {
    allInfo %<>%
      rbind(extra6)
  }
  extra6 = extra6 %>%
    select(sample_id, Run) %>%
    rename(v2_sample_id = sample_id) %>%
    inner_join(allIDs, by = join_by(Run == name), relationship = "many-to-many") %>%
    rename(v3_sample_id = sample_id, name = Run) %>%
    distinct(v2_sample_id, v3_sample_id, .keep_all = TRUE) %>%
    select(v2_sample_id, v3_sample_id, name)
  initialMatch %<>%
    rbind(extra6)
  mahaFile$status[mahaFile$Sample_ID %in% extra6$v2_sample_id] = "matched after OSR"

  ## THESE ONES SEEM TO BE LOST TO FOLLOW-UP CONFIRMATION
  Tab7a = read_csv("SraRunTableMOLDOVA.csv", show_col_types = FALSE, guess_max = 1e5) %>%
    select(Run, BioSample, Experiment, `Library Name`) %>%
    mutate_all(trimws) %>%
    mutate(Sample_Name = toupper(`Library Name`)) %>%
    select(-`Library Name`)
  extra7a = mahaUnmatched %>%
    inner_join(Tab7a, by = join_by(short_name == Sample_Name)) 
  mahaFile$status[mahaFile$Sample_ID %in% extra7a$sample_id] = "found in MOLDOVA"
  if (keepInfo) {
    allInfo %<>%
      rbind(extra7a)
  }
  extra7a = extra7a %>%
    select(sample_id, Run) %>%
    rename(v2_sample_id = sample_id) %>%
    inner_join(allIDs, by = join_by(Run == name), relationship = "many-to-many") %>%
    rename(v3_sample_id = sample_id, name = Run) %>%
    distinct(v2_sample_id, v3_sample_id, .keep_all = TRUE) %>%
    select(v2_sample_id, v3_sample_id, name)
  initialMatch %<>%
    rbind(extra7a)
  mahaFile$status[mahaFile$Sample_ID %in% extra7a$v2_sample_id] = "matched after MOLDOVA"
  
  ## THESE ONES SEEM TO BE ALREADY LOST FROM MAHA'S FILE
  Tab7b = read_csv("SraRunTableESTONIA.csv", show_col_types = FALSE, guess_max = 1e5) %>%
    select(Run, BioSample, Experiment, `Library Name`) %>%
    mutate_all(trimws) %>%
    mutate(Sample_Name = sapply(str_split(toupper(`Library Name`), "_"), dplyr::first)) %>%
    select(-`Library Name`)
  extra7b = mahaUnmatched %>%
    inner_join(Tab7b, by = join_by(short_name == Sample_Name)) 
  mahaFile$status[mahaFile$Sample_ID %in% extra7b$sample_id] = "found in ESTONIA"
  if (keepInfo) {
    allInfo %<>%
      rbind(extra7b)
  }
  extra7b = extra7b %>%
    select(sample_id, Run) %>%
    rename(v2_sample_id = sample_id) %>%
    inner_join(allIDs, by = join_by(Run == name), relationship = "many-to-many") %>%
    rename(v3_sample_id = sample_id, name = Run) %>%
    distinct(v2_sample_id, v3_sample_id, .keep_all = TRUE) %>%
    select(v2_sample_id, v3_sample_id, name)
  initialMatch %<>%
    rbind(extra7b)
  mahaFile$status[mahaFile$Sample_ID %in% extra7b$v2_sample_id] = "matched after ESTONIA"
  
  ## NOTE THAT THE FORMAT OF SEQ&TREAT IDENTFIERS IS TB0[4 DIGITS]0[3 DIGITS], AND THE MIDDLE 0 IS MISSING IN MAHA'S FILE!
  Tab8 = read_csv("SraRunTableSEQTREAT.csv", show_col_types = FALSE, guess_max = 1e5) %>%
    select(Run, BioSample, Experiment, `Library Name`) %>%
    mutate_all(trimws) %>%
    mutate(Sample_Name = paste0(str_sub(`Library Name`, 1, 7), str_sub(`Library Name`, 9, 11))) %>%
    select(-`Library Name`)
  extra8 = mahaUnmatched %>%
    inner_join(Tab8, by = join_by(short_name == Sample_Name)) 
  mahaFile$status[mahaFile$Sample_ID %in% extra8$sample_id] = "found in SEQTREAT"
  if (keepInfo) {
    allInfo %<>%
      rbind(extra8)
  }
  extra8 = extra8 %>%
    select(sample_id, Run) %>%
    rename(v2_sample_id = sample_id) %>%
    inner_join(allIDs, by = join_by(Run == name), relationship = "many-to-many") %>%
    rename(v3_sample_id = sample_id, name = Run) %>%
    distinct(v2_sample_id, v3_sample_id, .keep_all = TRUE) %>%
    select(v2_sample_id, v3_sample_id, name)
  initialMatch %<>%
    rbind(extra8)
  mahaFile$status[mahaFile$Sample_ID %in% extra8$v2_sample_id] = "matched after SEQTREAT"
  
  ## THESE ONES SEEM TO BE LOST TO FOLLOW-UP CONFIRMATION
  Tab9 = read_csv("SraRunTableSWISS.csv", show_col_types = FALSE, guess_max = 1e5) %>%
    select(Run, BioSample, Experiment, `Library Name`) %>%
    mutate_all(trimws) %>%
    mutate(Sample_Name = paste0(sapply(str_split(`Library Name`, ". "), dplyr::first), "_", sapply(str_split(`Library Name`, ". "), dplyr::nth, -2))) %>%
    select(-`Library Name`)
  extra9 = mahaUnmatched %>%
    inner_join(Tab9, by = join_by(short_name == Sample_Name))
  mahaFile$status[mahaFile$Sample_ID %in% extra9$sample_id] = "found in SWISS"
  if (keepInfo) {
    allInfo %<>%
      rbind(extra9)
  }
  extra9 = extra9 %>%
    select(sample_id, Run) %>%
    rename(v2_sample_id = sample_id) %>%
    inner_join(allIDs, by = join_by(Run == name), relationship = "many-to-many") %>%
    rename(v3_sample_id = sample_id, name = Run) %>%
    distinct(v2_sample_id, v3_sample_id, .keep_all = TRUE) %>%
    select(v2_sample_id, v3_sample_id, name)
  initialMatch %<>%
    rbind(extra9)
  mahaFile$status[mahaFile$Sample_ID %in% extra9$v2_sample_id] = "matched after SWISS"
  
  Tab9b = read_csv("SraRunTableARGENTINA.csv", show_col_types = FALSE, guess_max = 1e5) %>%
    select(Run, BioSample, Experiment, isolate) %>%
    mutate_all(trimws) %>%
    mutate(Sample_Name = isolate) %>%
    select(-isolate)
  extra9b = mahaUnmatched %>%
    inner_join(Tab9b, by = join_by(short_name == Sample_Name))
  mahaFile$status[mahaFile$Sample_ID %in% extra9b$sample_id] = "found in ARGENTINA"
  if (keepInfo) {
    allInfo %<>%
      rbind(extra9b)
  }
  extra9b = extra9b %>%
    select(sample_id, Run) %>%
    rename(v2_sample_id = sample_id) %>%
    inner_join(allIDs, by = join_by(Run == name), relationship = "many-to-many") %>%
    rename(v3_sample_id = sample_id, name = Run) %>%
    distinct(v2_sample_id, v3_sample_id, .keep_all = TRUE) %>%
    select(v2_sample_id, v3_sample_id, name)
  initialMatch %<>%
    rbind(extra9b)
  mahaFile$status[mahaFile$Sample_ID %in% extra9b$v2_sample_id] = "matched after ARGENTINA"
  
  # bTab1 = read_csv("SraRunTableRUSSIA.csv") %>%
  #   select(BioSample, Run, Experiment)
  # b1 = mahaUnmatchedBio %>%
  #   inner_join(bTab1, by = join_by(name == BioSample))
  # mahaFile$status[mahaFile$Sample_ID %in% b1$sample_id] = "found by BioSample in RUSSIA"
  # b1 = b1 %>%
  #   select(sample_id, Run, Experiment) %>%
  #   rename(v2_sample_id = sample_id) %>%
  #   inner_join(allIDs, by = join_by(Run == name), relationship = "many-to-many") %>%
  #   rename(v3_sample_id = sample_id, name = Run) %>%
  #   distinct(v2_sample_id, v3_sample_id, .keep_all = TRUE) %>%
  #   select(v2_sample_id, v3_sample_id, name)
  # mahaFile$status[mahaFile$Sample_ID %in% b1$v2_sample_id] = "matched by BioSample in RUSSIA"
  # 
  # bTab2 = read_csv("SraRunTableOCICB.csv") %>%
  #   select(BioSample, Run, Experiment)
  # b2 = mahaUnmatchedBio %>%
  #   inner_join(bTab2, by = join_by(name == BioSample))
  # mahaFile$status[mahaFile$Sample_ID %in% b2$sample_id] = "found by BioSample in OCICB"
  # b2 = b2 %>%
  #   select(sample_id, Run, Experiment) %>%
  #   rename(v2_sample_id = sample_id) %>%
  #   inner_join(allIDs, by = join_by(Run == name), relationship = "many-to-many") %>%
  #   rename(v3_sample_id = sample_id, name = Run) %>%
  #   distinct(v2_sample_id, v3_sample_id, .keep_all = TRUE) %>%
  #   select(v2_sample_id, v3_sample_id, name)
  # mahaFile$status[mahaFile$Sample_ID %in% b2$v2_sample_id] = "matched by BioSample in OCICB"
  # 
  # bTab3 = read_csv("SraRunTableSANGER.csv") %>%
  #   select(BioSample, Run, Experiment)
  # b3 = mahaUnmatchedBio %>%
  #   inner_join(bTab3, by = join_by(name == BioSample))
  # mahaFile$status[mahaFile$Sample_ID %in% b3$sample_id] = "found by BioSample in SANGER"
  # b3 = b3 %>%
  #   select(sample_id, Run, Experiment) %>%
  #   rename(v2_sample_id = sample_id) %>%
  #   inner_join(allIDs, by = join_by(Run == name), relationship = "many-to-many") %>%
  #   rename(v3_sample_id = sample_id, name = Run) %>%
  #   distinct(v2_sample_id, v3_sample_id, .keep_all = TRUE) %>%
  #   select(v2_sample_id, v3_sample_id, name)
  # mahaFile$status[mahaFile$Sample_ID %in% b3$v2_sample_id] = "matched by BioSample in SANGER"
  # 
  # bTab4 = read_csv("SraRunTableSERBIA.csv") %>%
  #   select(BioSample, Run, Experiment)
  # b4 = mahaUnmatchedBio %>%
  #   inner_join(bTab4, by = join_by(name == BioSample))
  # mahaFile$status[mahaFile$Sample_ID %in% b4$sample_id] = "found by BioSample in SERBIA"
  # b4 = b4 %>%
  #   select(sample_id, Run, Experiment) %>%
  #   rename(v2_sample_id = sample_id) %>%
  #   inner_join(allIDs, by = join_by(Run == name), relationship = "many-to-many") %>%
  #   rename(v3_sample_id = sample_id, name = Run) %>%
  #   distinct(v2_sample_id, v3_sample_id, .keep_all = TRUE) %>%
  #   select(v2_sample_id, v3_sample_id, name)
  # mahaFile$status[mahaFile$Sample_ID %in% b4$v2_sample_id] = "matched by BioSample in SERBIA"
  
  # leftOver = mahaUnmatchedBio %>% 
  #   filter(sample_id %in% (mahaFile %>% filter(status == "unmatched") %>% pull(Sample_ID))) %>% 
  #   pull(name) %>%
  #   sort()
  # res = vector("list", length(leftOver))
  # for (ind in 1:length(leftOver)) { 
  #   res[[ind]] = retrieve_sra_metadata(leftOver[ind])
  #   print(ind)
  # }
  # names(res) = leftOver
  
  ### TODO: CONSIDER SKIPPING THIS IF keepInfo IS FALSE!
  results = lookupDriver(inputTab = mahaUnmatchedBio) %>%
    as_tibble()
  newResults = results %>%
    inner_join(mahaUnmatchedBio, by = join_by(sample_accession == name)) 
  mahaFile$status[mahaFile$BioSample %in% newResults$sample_accession] = "found by BioSample"
  newResults = newResults %>%
    inner_join(allIDs, by = join_by(run_accession == name), relationship = "many-to-many") %>%
    rename(v2_sample_id = sample_id.x, v3_sample_id = sample_id.y, name = run_accession) %>%
    distinct(v2_sample_id, v3_sample_id, .keep_all = TRUE) %>%
    select(v2_sample_id, v3_sample_id, name)
  initialMatch %<>%
    rbind(newResults)
  mahaFile$status[mahaFile$Sample_ID %in% newResults$v2_sample_id] = "matched by BioSample"
  
  extraChecks = results %>%
    filter(is.na(run_accession)) %>%
    pull(sample_accession)
  
  write_csv(initialMatch, file = "InitialMatch.csv")
  if (keepInfo) {
    mahaFileExtra = mahaFile %>% 
      inner_join(allInfo, by = join_by(Sample_ID == sample_id))
    mahaFileExtra = mahaFileExtra %>%
      select(-BioSample.x) %>%
      rename(BioSample = BioSample.y)
    write_csv(mahaFileExtra, file = "ExtraInfoV2.csv")
  }
  # problemFile = mahaIDs %>%
  #   filter(!(name %in% allIDs$name)) %>%
  #   group_by(sample_id) %>%
  #   mutate(N  = n()) %>%
  #   filter(N == 2) %>%
  #   ungroup()
  # missingSampleIDs = missingIDs[str_detect(missingIDs, ACCESSIONS[["BIOSAMPLE"]])]
  # missingOtherIDs  = base::setdiff(missingIDs, missingSampleIDs)
  setwd(initDir)
  initialMatch
}

## For the quality score analysis:

## This should produce more or less identical results to the ones we are seeing now

## A0=Sys.time()
## A1K = mainDriver(EXTRACTION_ID = EXTRA_ID, minMAF = MAF_THRESHOLD_REGULAR, minQ = QUALITY_THRESHOLD_STRICT, fast = TRUE)
## B0=Sys.time()-A0

## Use this code to compute regular and relaxed variant catalogs with different quality score thresholds:

# Results = vector("list")
# for (threshold in c(250, 500, 1000)) {
# 	for (relaxed in c(FALSE, TRUE)) {
# 		curName = paste0(ifelse(relaxed, "Relaxed", "Regular"), "Q", threshold)	
# 		print(curName)
# 		Results[[curName]] = mainDriver(EXTRACTION_ID = EXTRA_ID, 
# 		minMAF = ifelse(relaxed, MAF_THRESHOLD_RELAXED, MAF_THRESHOLD_REGULAR), minQ = threshold, fast = TRUE) 
# 	}
# }

## Use this code to compute regular and relaxed results (sensitivity and specificity) with different quality score thresholds:

# Results = vector("list")
# for (threshold in c(250, 500, 1000)) {
# 	for (relaxed in c(FALSE, TRUE)) {
# 		curName = paste0(ifelse(relaxed, "Relaxred", "Regular"), "Q", threshold)	
# 		Results[[curName]] = computeSensSpec(version = 2, relaxed = relaxed, minQ = threshold, EXTRACTION_ID = EXTRA_ID, useOrphanData = FALSE)
# 	}
# }

## Use this code to compute extended information about compensatory mutations:

## W = computeSensSpec(skipCompensatory = FALSE, recent = TRUE, useOrphanData = FALSE, EXTRACTION_ID = EXTRA_ID)

# Compare two FASTQ files ignoring identifiers and read order (generated by Claude)
# Separately compares sequences and qualities to distinguish different scenarios
# Usage example: result <- compare_fastq_files("file1.fastq", "file2.fastq")
compare_fastq_files <- function(file1, file2) {
  # Read both files
  lines1 <- readLines(file1)
  lines2 <- readLines(file2)
  
  # Check if files have same number of lines
  if (length(lines1) != length(lines2)) {
    cat("Files have different number of lines\n")
    return(list(sequences_identical = FALSE, qualities_identical = FALSE))
  }
  
  # FASTQ format: every 4 lines = 1 read (identifier, sequence, +, quality)
  n_reads <- length(lines1) / 4
  
  if (length(lines1) %% 4 != 0) {
    cat("Invalid FASTQ format - line count not divisible by 4\n")
    return(list(sequences_identical = FALSE, qualities_identical = FALSE))
  }
  
  # Extract sequences and qualities separately from both files
  sequences1 <- character(n_reads)
  qualities1 <- character(n_reads)
  sequences2 <- character(n_reads)
  qualities2 <- character(n_reads)
  
  for (i in 1:n_reads) {
    seq_line <- 4 * (i - 1) + 2  # sequence line
    qual_line <- 4 * (i - 1) + 4  # quality line
    
    sequences1[i] <- lines1[seq_line]
    qualities1[i] <- lines1[qual_line]
    sequences2[i] <- lines2[seq_line]
    qualities2[i] <- lines2[qual_line]
  }
  
  # Sort sequences and qualities separately for comparison
  sequences1_order  <- order(sequences1)
  sequences1_sorted <- sequences1[sequences1_order]
  qualities1_sorted <- qualities1[sequences1_order]
  
  sequences2_order  <- order(sequences2)
  sequences2_sorted <- sequences2[sequences2_order]
  qualities2_sorted <- qualities2[sequences2_order]
  
  # Compare sequences and qualities independently
  sequences_identical <- identical(sequences1_sorted, sequences2_sorted)
  qualities_identical <- identical(qualities1_sorted, qualities2_sorted)
  
  # Report results
  cat("Total reads:", n_reads, "\n")
  cat("Sequences identical (ignoring order):", sequences_identical, "\n")
  cat("Qualities identical (ignoring order):", qualities_identical, "\n")
  
  # Determine scenario
  if (sequences_identical && qualities_identical) {
    cat("Result: Files are identical (same sequences with same qualities)\n")
  } else if (sequences_identical && !qualities_identical) {
    cat("Result: Sequences are identical, but qualities differ\n")
  } else if (!sequences_identical && qualities_identical) {
    cat("Result: Qualities are identical, but sequences differ\n")
  } else {
    cat("Result: Both sequences and qualities differ\n")
  }
  
  return(list(sequences_identical = sequences_identical, 
              qualities_identical = qualities_identical))
}

# Biosample ID Lookup Script (generated by Claude)
# Retrieves runs and experiments for ENA and SRA biosample IDs

# Function to determine if ID is ENA or SRA format
identify_database <- function(biosample_id) {
  if (grepl("^SAMN", biosample_id)) {
    return("SRA")
  } else if (grepl("^SAMEA", biosample_id)) {
    return("ENA")
  } else if (grepl("^SAME", biosample_id)) {
    return("ENA")  # ENA also uses SAME prefix
  } else {
    return("Unknown")
  }
}

# Function to query ENA API for biosample information
query_ena_biosample <- function(biosample_id) {
  base_url <- "https://www.ebi.ac.uk/ena/portal/api/search"
  
  # Query for experiments linked to this biosample
  exp_query <- paste0(base_url, 
                      "?result=read_run",
                      "&query=sample_accession=", 
                      biosample_id,
                      "&format=json",
                      "&fields=experiment_accession,run_accession,sample_accession,study_accession,",
                      "fastq_ftp,fastq_md5,submitted_ftp,submitted_md5,sra_ftp,sra_md5")
  
  tryCatch({
    response <- GET(exp_query)
    if (status_code(response) == 200) {
      content <- content(response, as = "text", encoding = "UTF-8")
      if (content != "") {
        data <- fromJSON(content)
        if (nrow(data) > 0) {
          return(data)
        }
      }
    }
    return(data.frame())
  }, error = function(e) {
    message(paste("Error querying ENA for", biosample_id, ":", e$message))
    return(data.frame())
  })
}

# Function to query SRA via NCBI E-utilities
query_sra_biosample <- function(biosample_id) {
  # First, search for SRA experiments linked to this biosample
  search_url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
                       "?db=sra",
                       "&term=", biosample_id, "[BioSample]",
                       "&retmode=json")
  
  tryCatch({
    search_response <- GET(search_url)
    if (status_code(search_response) == 200) {
      search_data <- fromJSON(content(search_response, as = "text"))
      
      if (length(search_data$esearchresult$idlist) > 0) {
        # Get detailed information for found IDs
        ids <- paste(search_data$esearchresult$idlist, collapse = ",")
        
        fetch_url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                            "?db=sra",
                            "&id=", ids,
                            "&rettype=xml")
        Sys.sleep(0.25)  # To respect rate limits
        fetch_response <- GET(fetch_url)
        if (status_code(fetch_response) == 200) {
          xml_content <- read_xml(content(fetch_response, as = "text"))
          
          # Parse XML to extract experiment and run information
          experiments <- xml_find_all(xml_content, ".//EXPERIMENT")
          
          results <- map_dfr(experiments, function(exp) {
            exp_acc <- xml_attr(exp, "accession")
            
            # Find associated runs
            exp_package <- xml_parent(xml_parent(exp))
            runs <- xml_find_all(exp_package, ".//RUN")
            
            # Get study accession
            study <- xml_find_first(exp_package, ".//STUDY")
            study_acc <- xml_attr(study, "accession")
            if (length(runs) > 0) {
              run_data <- map_dfr(runs, function(run) {
                run_acc <- xml_attr(run, "accession")
                # Extract file information including MD5 checksums
                sra_files <- xml_find_all(run, ".//SRAFile")
                
                if (length(sra_files) > 0) {
                  file_info <- map_dfr(sra_files, function(file) {
                    data.frame(
                      experiment_accession = exp_acc,
                      run_accession = run_acc,
                      sample_accession = biosample_id,
                      study_accession = study_acc,
                      filename = xml_attr(file, "filename"),
                      file_url = xml_attr(file, "url"),
                      file_md5 = xml_attr(file, "md5"),
                      file_size = xml_attr(file, "size"),
                      stringsAsFactors = FALSE
                    )
                  })
                  return(file_info)
                } else {
                  # No file information available
                  data.frame(
                    experiment_accession = exp_acc,
                    run_accession = run_acc,
                    sample_accession = biosample_id,
                    study_accession = study_acc,
                    filename = NA,
                    file_url = NA,
                    file_md5 = NA,
                    file_size = NA,
                    stringsAsFactors = FALSE
                  )
                }
              })
              return(run_data)
            } else {
              data.frame(
                experiment_accession = exp_acc,
                run_accession = NA,
                sample_accession = biosample_id,
                study_accession = NA,
                filename = NA,
                file_url = NA,
                file_md5 = NA,
                file_size = NA,
                stringsAsFactors = FALSE
              )
            }
          })
          return(results)
        }
      }
    }
    return(data.frame())
  }, error = function(e) {
    message(paste("Error querying SRA for", biosample_id, ":", e$message))
    return(data.frame())
  })
}

# Main function to process a list of biosample IDs
process_biosample_ids <- function(biosample_ids) {
  numIDs <- length(biosample_ids)
  print(paste("There are", numIDs, "biosample IDs to process."))
  results <- map_dfr(1:numIDs, function(index) {
    if (index %% 10 == 0) {
      message(paste("Processing ID number", index))
    }
    Sys.sleep(0.25)  # To respect rate limits
    id <- biosample_ids[index]
    db_type <- identify_database(id)
    if (db_type == "ENA") {
      data <- query_ena_biosample(id)
    } else if (db_type == "SRA") {
      data <- query_sra_biosample(id)
    } else {
      message(paste("Unknown biosample format for", id, ". Trying both databases..."))
      # Try both databases
      ena_data <- query_ena_biosample(id)
      sra_data <- query_sra_biosample(id)
      data <- rbind(ena_data, sra_data)
    }
    
    if (nrow(data) == 0) {
      # Return a row with NA values to indicate no results found
      data.frame(
        experiment_accession = NA,
        run_accession = NA,
        sample_accession = id,
        study_accession = NA,
        fastq_ftp = NA,
        fastq_md5 = NA,
        submitted_ftp = NA,
        submitted_md5 = NA,
        sra_ftp = NA,
        sra_md5 = NA,
        filename = NA,
        file_url = NA,
        file_md5 = NA,
        file_size = NA,
        database = db_type,
        stringsAsFactors = FALSE
      )
    } else {
      data$database <- db_type
      return(data)
    }
  })
  
  return(results)
}

lookupDriver = function(inputTab) {
  # Define your list of biosample IDs
  biosample_ids <- inputTab$name
  
  # Process the biosample IDs
  message("Starting biosample lookup...")
  results <- process_biosample_ids(biosample_ids)
  
  # Display results
  print(results)
  
  # Save results to CSV
  write.csv(results, "biosample_lookup_results.csv", row.names = FALSE)
  message("Results saved to biosample_lookup_results.csv")
  
  # Summary statistics
  cat("\nSummary:\n")
  cat("Total biosamples processed:", length(unique(results$sample_accession)), "\n")
  cat("Total experiments found:", sum(!is.na(results$experiment_accession)), "\n")
  cat("Total runs found:", sum(!is.na(results$run_accession)), "\n")
  
  # Show samples with no results
  no_results <- results[is.na(results$experiment_accession), "sample_accession"]
  if (length(no_results) > 0) {
    cat("Biosamples with no results found:", paste(no_results, collapse = ", "), "\n")
  }
  results
}

getMD5Sums = function() {
  Tab = read_csv("~/Downloads/md5sum.csv") %>% ### Removing the other prep strategies as they don't enter v3
    filter(!is.na(sample_id) & algorithm == "MD5" & library_preparation_strategy %in% c("WGA", "WGS")) %>%
    select(-id_y) ### Removing these because they are identical to sequencing_id
  Tab = Tab %>%
    select(value, sequencing_data_id, filename, file_size, library_name, sample_id) %>%
    rename(MD5 = value, Seq_ID = sequencing_data_id, Filename = filename, File_Size = file_size, Lib_Name = library_name, Sample_ID = sample_id)
  
}

matchBiosamples = function() {
  mahaFile = read_xlsx("~/Downloads/FINDConsulting/Maha Regression Supplementary File.xlsx") %>%
    select(Sample_ID, BioSample, Sample_Name, Lineage) %>%
    mutate_all(trimws) %>%
    mutate_at("Sample_ID", as.numeric)
  ### Preprocess by splitting compound IDs into just the SRR/ERR part
  mahaFileProcessed = mahaFile %>%
    filter(!is.na(BioSample)) %>%
    mutate_at("BioSample", ~ifelse(str_detect(., "_"), str_split_fixed(., "_", n = 2)[,2], .))
  infoTab = read_csv("~/Downloads/FINDConsulting/DatabaseElements/submission_samplealias.csv") %>%
    select(-verdicts)
  allMatches = inner_join(mahaFileProcessed, infoTab, by = join_by(BioSample == name)) %>%
    rename(v2_sample_id = Sample_ID, v3_sample_id = sample_id)
  reducedMatches = allMatches %>%
    select(v2_sample_id, v3_sample_id, BioSample) %>%
    distinct(v2_sample_id, v3_sample_id, .keep_all = TRUE)
  otherMatches = read_csv("~/Downloads/FINDConsulting/SOLO Algorithm Input Files/DATABASE EXTRACTION files/2025-11-04T13-45-18.233070_jr_c34fade82bec882a90024c6d0a96405749a52aaeab3286803e85a6ef9ee92d50/identifiers/InitialMatch.csv") 
  ### Check consistency between the two files
  check1 = reducedMatches %>%
    inner_join(otherMatches, by = join_by(v2_sample_id)) %>%
    filter(v3_sample_id.x != v3_sample_id.y)
  stopifnot(nrow(check1) == 0)
  check2 = reducedMatches %>%
    inner_join(otherMatches, by = join_by(v3_sample_id), relationship = "many-to-many") %>%
    group_by(v3_sample_id) %>%
    mutate(N1 = n_distinct(v2_sample_id.x), N2 = n_distinct(v2_sample_id.y)) %>%
    filter(N1 != N2 & v2_sample_id.x != v2_sample_id.y)
  check3 = reducedMatches %>% 
    group_by(v3_sample_id) %>% 
    mutate(N = n_distinct(v2_sample_id)) %>% 
    filter(N > 1) %>% 
    ungroup() %>% 
    arrange(v3_sample_id)
  ### Preprocess the remaining problem IDs by splitting compound IDs into just the BioSample part
  mahaFileOther = mahaFile %>%
    filter(!is.na(BioSample) & !(Sample_ID %in% reducedMatches$v2_sample_id)) %>%
    mutate(BioSample = str_split_fixed(Sample_Name, "_", n = 2)[,1]) %>%
    inner_join(infoTab, by = join_by(BioSample == name)) %>%
    rename(v2_sample_id = Sample_ID, v3_sample_id = sample_id) %>%
    select(v2_sample_id, v3_sample_id, BioSample)
  jodyFile = read_csv("~/Downloads/FINDConsulting/jody/sample_match_attempt.csv") %>%
    select(Sample_ID, New_ID, New_ID_BioSample) %>%
    inner_join(infoTab, by = join_by(New_ID_BioSample == name)) %>%
    rename(v2_sample_id = Sample_ID, v3_sample_id = sample_id, BioSample = New_ID_BioSample) %>%
    select(v2_sample_id, v3_sample_id, BioSample)
  otherJodyFile = read_csv("~/Downloads/FINDConsulting/jody/additional_matching_attempt.csv") %>%
    select(Sample_ID, v3_sample_id, BioSample) %>%
    inner_join(infoTab, by = join_by(BioSample == name)) %>%
    rename(v2_sample_id = Sample_ID) %>%
    select(v2_sample_id, v3_sample_id, BioSample)
  otherMatches = otherMatches %>%
    rename(BioSample = name) %>%
    filter(str_detect(BioSample, ACCESSIONS[["BIOSAMPLE"]]))
  output = reducedMatches %>%
    bind_rows(otherMatches) %>%
    bind_rows(mahaFileOther) %>%
    bind_rows(jodyFile) %>%
    bind_rows(otherJodyFile) %>%
    distinct(v2_sample_id, v3_sample_id, .keep_all = TRUE) %>%
    arrange(v2_sample_id)
  extraData = read_csv("~/Downloads/FINDConsulting/jody/SRATables/ExtraInfoV2.csv")
  extraOutput = extraData %>%
    filter(!(Sample_ID %in% output$v2_sample_id)) %>%
    inner_join(infoTab, by = join_by(BioSample == name), relationship = "many-to-many") %>%
    select(-name) %>%
    rename(v2_sample_id = Sample_ID, v3_sample_id = sample_id) %>%
    select(v2_sample_id, v3_sample_id, BioSample)
  output = output %>% 
    bind_rows(extraOutput) %>%
    distinct(v2_sample_id, v3_sample_id, .keep_all = TRUE) %>%
    arrange(v2_sample_id)
  write_csv(output, file = "SOLO Algorithm Input Files/DATABASE EXTRACTION files/2025-11-04T13-45-18.233070_jr_c34fade82bec882a90024c6d0a96405749a52aaeab3286803e85a6ef9ee92d50/identifiers/ExtendedMatchBioSample.csv")
}

### Check about reconsenting issues next?!

prepMatchTab = function() {
  matchTab = read_csv("jody/ExtendedMatchBioSample.csv",
                      guess_max = LARGE_NUMBER, show_col_types = FALSE)
  badTab = read_csv("jody/DodgyMatches.csv")
  matchTab = matchTab %>% 
    anti_join(badTab)
  matchTab = matchTab %>%
    group_by(v2_sample_id) %>%
    mutate(N = n()) %>%
    filter(N == 1) %>%
    ungroup() %>%
    select(-N)
  matchTab = matchTab %>%
    group_by(v3_sample_id) %>%
    mutate(N = n()) %>%
    filter(N == 1) %>%
    ungroup() %>%
    select(-N)
  matchTab
}

tableForJody = function() {
  matchTab = prepMatchTab()
  v2Tab = read_csv("mutation-catalogue-2023/Results/2023-04-25T06_00_10.443990_jr_b741dc136e079fa8583604a4915c0dc751724ae9880f06e7c2eacc939e086536/CompleteDataset.csv", 
                   guess_max = LARGE_NUMBER, show_col_types = FALSE) %>%
    select(sample_id, drug, phenotype) %>%
    distinct(sample_id, drug, .keep_all = TRUE) %>%
    rename(v2_sample_id = sample_id)
  v2Tab = v2Tab %>%
    inner_join(matchTab, by = join_by(v2_sample_id)) %>%
    arrange(v3_sample_id, drug)
  v2Tab = v2Tab %>%
    rename(`Sample ID` = BioSample)
  v2Tab = v2Tab %>%
    select(-v2_sample_id, -v3_sample_id)
  v2Tab = v2Tab %>%
    mutate(`FASTQ prefix` = NA, Country = NA, `Sampling date` = NA, `DST Method` = NA)
  v2Tab = v2Tab %>%
    mutate(drug = SHORT_NAMES[match(drug, DRUG_LIST)])
  v2TabWide = v2Tab %>%
    pivot_wider(names_from = drug, values_from = phenotype)
  write_excel_csv(v2TabWide, file = "table_for_jody.xlsx")
  write_csv(v2TabWide, file = "table_for_jody.csv")
  v2Tab
}

readV2Files = function(extra = TRUE) {
  drugTab = read_parquet('V2Archive/genphensql.drug/part-00000-20af102a-eec2-4d63-8b96-f4043a7747fa-c000.gz.parquet')
  mediumTab = read_parquet('V2Archive/genphensql.growth_medium/part-00000-9342a4be-4c6c-40bc-afec-066393707cb0-c000.gz.parquet')
  methodTab = read_parquet('V2Archive/genphensql.phenotypic_drug_susceptibility_assessment_method/part-00000-9a8b7222-9b26-4332-86d9-f28d290dedab-c000.gz.parquet')
  resultTab = read_parquet('V2Archive/genphensql.phenotypic_drug_susceptibility_test/part-00000-e1894afc-6f6e-4788-8a07-ee806628bae3-c000.gz.parquet')
  resultTab = resultTab %>%
    left_join(drugTab) %>%
    left_join(mediumTab) %>%
    left_join(methodTab)
  resultTab = resultTab %>%
    select(sample_id, drug_name, medium_name, method_name, concentration, test_result)
  matchTab = prepMatchTab()
  write_csv(matchTab, file = "jody/FinalMatches.csv")
  resultTab = resultTab %>%
    rename(v2_sample_id = sample_id) %>%
    inner_join(matchTab)
  ####### NEWLY ADDED NOV 26, 2025 #########
  if (extra) {
    extraResultTab = resultTab %>%
      mutate(drug = EXTRA_NAMES[match(drug_name, EXTRA_LIST)]) %>%
      select(-drug_name) %>%
      filter(!is.na(drug)) %>%
      mutate(drug_CC = paste0(drug, ifelse(is.na(concentration), '', paste0(' (', concentration, ')')))) %>%
      select(-drug, -concentration) %>%
      distinct() %>%
      mutate(`FASTQ prefix` = NA, Country = NA, `Sampling date` = NA) %>%
      rename(`DST Method` = medium_name, `Assessment method` = method_name, `Sample ID` = BioSample) %>%
      select(-v2_sample_id, -v3_sample_id) %>%
      select(`Sample ID`, `FASTQ prefix`, Country, `Sampling date`, `DST Method`, `Assessment method`, everything()) %>%
      distinct()
    extraWideTab = extraResultTab %>%
      pivot_wider(names_from = drug_CC, values_from = test_result)
    extraSplitTab = extraWideTab %>%
      split(extraWideTab$`DST Method`)
    for (ind in 1:length(extraSplitTab)) {
      curName = names(extraSplitTab)[ind]
      write_csv(extraSplitTab[[ind]], file = paste0('jody/Extra_DST_table_v2_', curName, '.csv'))
    }
  }
  #######
  resultTab = resultTab %>%
    mutate(drug = SHORT_NAMES[match(drug_name, DRUG_LIST)]) %>%
    select(-drug_name) %>%
    filter(!is.na(drug))
  resultTab = resultTab %>%
    mutate(drug_CC = paste0(drug, ifelse(is.na(concentration), '', paste0(' (', concentration, ')'))))
  resultTab = resultTab %>%
    select(-drug, -concentration)
  resultTab = distinct(resultTab)
  intermedTab = resultTab %>%
    summarise(n = dplyr::n(), .by = c(v2_sample_id, medium_name, method_name, v3_sample_id, BioSample, drug_CC))
  contradictoryResults = intermedTab %>%
    filter(n > 1)
  save(contradictoryResults, file = "jody/ContradictoryDST.csv")
  resultTab = resultTab %>%
    anti_join(contradictoryResults, by = join_by(BioSample, drug_CC, medium_name))
  resultTab = resultTab %>%
    mutate(`FASTQ prefix` = NA, Country = NA, `Sampling date` = NA) %>%
    rename(`DST Method` = medium_name, `Assessment method` = method_name, `Sample ID` = BioSample) %>%
    select(-v2_sample_id, -v3_sample_id) %>%
    select(`Sample ID`, `FASTQ prefix`, Country, `Sampling date`, `DST Method`, `Assessment method`, everything()) %>%
    distinct()
  wideTab = resultTab %>%
    pivot_wider(names_from = drug_CC, values_from = test_result)
  splitTab = wideTab %>%
    split(wideTab$`DST Method`)
  for (ind in 1:length(splitTab)) {
    curName = names(splitTab)[ind]
    write_csv(splitTab[[ind]], file = paste0('jody/DST_table_v2_', curName, '.csv'))
  }
  splitTab
}


readV2MICFiles = function(extra = TRUE) {
  drugTab = read_parquet('V2Archive/genphensql.drug/part-00000-20af102a-eec2-4d63-8b96-f4043a7747fa-c000.gz.parquet')
  resultTab = read_parquet('V2Archive/genphensql.minimum_inhibitory_concentration_test/part-00000-e2397344-c9e0-42d6-a386-eee778e2514b-c000.gz.parquet')
  resultTab = resultTab %>%
    left_join(drugTab) 
  resultTab = resultTab %>%
    select(sample_id, plate, drug_name, mic_value)
  # matchTab = prepMatchTab()
  # save(matchTab, file = "jody/FinalMatches.csv")
  load("jody/FinalMatches.csv")
  resultTab = resultTab %>%
    rename(v2_sample_id = sample_id) %>%
    inner_join(matchTab)
  ####### NEWLY ADDED NOV 26, 2025 #########
  if (extra) {
    extraResultTab = resultTab %>%
      mutate(drug = EXTRA_NAMES[match(drug_name, EXTRA_LIST)]) %>%
      select(-drug_name) %>%
      filter(!is.na(drug)) %>%
      distinct() %>%
      mutate(`FASTQ prefix` = NA) %>%
      rename(`DST Method` = plate, `Sample ID` = BioSample) %>%
      select(-v2_sample_id, -v3_sample_id) %>%
      select(`Sample ID`, `DST Method`, `FASTQ prefix`, everything()) %>%
      distinct()
    extraWideTab = extraResultTab %>%
      pivot_wider(names_from = drug, values_from = mic_value)
    extraSplitTab = extraWideTab %>%
      split(extraWideTab$`DST Method`)
    for (ind in 1:length(extraSplitTab)) {
      curName = names(extraSplitTab)[ind]
      write_csv(extraSplitTab[[ind]], file = paste0('jody/Extra_MIC_table_v2_', curName, '.csv'))
    }
  }
  #######
  resultTab = resultTab %>%
    mutate(drug = SHORT_NAMES[match(drug_name, DRUG_LIST)]) %>%
    select(-drug_name) %>%
    filter(!is.na(drug))
  
  # Don't need this I think
  ########
  # resultTab = resultTab %>%
  #   mutate(drug_CC = paste0(drug, ifelse(is.na(concentration), '', paste0(' (', concentration, ')'))))
  # resultTab = resultTab %>%
  #   select(-drug, -concentration)
  #############
  
  resultTab = distinct(resultTab)
  intermedTab = resultTab %>%
    summarise(n = dplyr::n(), .by = c(v2_sample_id, plate, v3_sample_id, BioSample, drug))
  contradictoryResults = intermedTab %>%
    filter(n > 1)
  save(contradictoryResults, file = "jody/ContradictoryMIC.csv")
  resultTab = resultTab %>%
    anti_join(contradictoryResults, by = join_by(BioSample, plate))
  resultTab = resultTab %>%
    mutate(`FASTQ prefix` = NA) %>%
    rename(`DST Method` = plate, `Sample ID` = BioSample) %>%
    select(-v2_sample_id, -v3_sample_id) %>%
    select(`Sample ID`, `DST Method`, `FASTQ prefix`, everything()) %>%
    distinct()
  wideTab = resultTab %>%
    pivot_wider(names_from = drug, values_from = mic_value)
  splitTab = wideTab %>%
    split(wideTab$`DST Method`)
  for (ind in 1:length(splitTab)) {
    curName = names(splitTab)[ind]
    write_csv(splitTab[[ind]], file = paste0('jody/MIC_table_v2_', curName, '.csv'))
  }
  splitTab
}
