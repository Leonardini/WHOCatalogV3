#' Execute the original SOLO algorithm (Walker et al., PMID 26116186)
#'
#' inputTab requires columns sample_id, mutation, phenotype and either SOnly or het.
#' maxIter is a positive integer or Inf. If removeSOnly = FALSE, the first stage is skipped.
#' If listIsolates = TRUE, all isolates classified as solo at each iteration are recorded.
#' @noRd
originalSOLO = function(inputTab, maxIter = MAX_ITER, removeSOnly = TRUE, listIsolates = FALSE) {
  useHet = ("het" %in% colnames(inputTab))
  ## Construct bipartite directed graph with each edge having a sample as tail and a mutation it contains as head
  G = graph_from_data_frame(inputTab, directed = TRUE)
  ## The mutations are in column 2 
  mutations = na.omit(unique(inputTab[[2]]))
  ## The output matrix contains, for each mutation, its classification, iteration classified, SOLO R/total counts
  output = matrix(NA, length(mutations), 4, dimnames = list(mutations, c("class", "iter", "Tcnt", "Rcnt")))
  iter = 1L
  currentSetS = c()
  if (listIsolates) {
    isoList = tibble()
  }
  ## The loop proceeds till either maxIter is reached, all mutations are classified or no new ones are classified S
  while (iter <= maxIter && length(mutations) > 0 && (iter == 1L || length(currentSetS) > 0)) {
    ## Remove mutations classified S from the list of active mutations and the graph, along with incident edges
    G %<>% delete_vertices(currentSetS)
    mutations %<>% setdiff(currentSetS)
    ## Compute the total out-degree, which is 0 for mutations and is the number of mutations contained for samples 
    Tdegree = igraph::degree(G, mode = "out")
    ## The singleton samples are those containing a single mutation
    currentSingletons = names(Tdegree)[Tdegree == 1]
    ## The initial set of mutations classified as S is empty
    currentSetS = c()
    ## Start by classifying all the mutations occurring only in S isolates as S at the first iteration, if necessary
    if (iter == 1L && removeSOnly) {
      if (useHet) {
        ## Use what we know; extract the subgraph containing only the R samples and non-het mutations
        Rsubgraph = igraph::subgraph_from_edges(G, eids = E(G)[phenotype == "R" & !het], delete.vertices = FALSE)
        ## Compute the degrees in this subgraph; degree 0 means a mutation never occurs in R samples as a non-het
        Rdegree = igraph::degree(Rsubgraph)
        ## Also extract the subgraph containing only the S samples and non-het mutations
        Ssubgraph = igraph::subgraph_from_edges(G, eids = E(G)[phenotype == "S" & !het], delete.vertices = FALSE)
        ## Compute the degrees in this subgraph; degree > 0 means a mutation occurs in some S samples as a non-het
        Sdegree = igraph::degree(Ssubgraph)
        ## Classify all the mutations satisfying both criteria as S
        currentSetS = mutations[Rdegree[mutations] == 0 & Sdegree[mutations] > 0]
      } else { 
        ## Extract the mutations marked SOnly; some are excluded due to hets
        currentSetS = inputTab %>% 
          dplyr::filter(SOnly) %>%
          pull(2) %>%
          unique()
      }
    }
    ## Compute the edges connecting the singleton samples to their mutations, and get their endpoints as a matrix
    if (useHet) {
      currentEdgeSeq = E(G)[.from(currentSingletons) & !het]
    } else {
      currentEdgeSeq = E(G)[.from(currentSingletons)]
    }
    currentEdges   = igraph::ends(G, currentEdgeSeq)
    ## The full counts are obtained by looking at all the mutations occurring in singletons (i.e. solo mutations)
    currentTCounts  = table(currentEdges[,2])
    currentSet      = names(currentTCounts)
    ## Identify the R edges among them (they originate from R samples) and classify their mutation endpoints as R
    currentREdges  = (currentEdgeSeq$phenotype == "R")
    currentSetR    = currentEdges[currentREdges, 2]
    ## The mutations found in R edges are counted as SOLO_R; note that after tabulation, the names become unique
    currentRCounts  = table(currentSetR)
    currentSetR     = names(currentRCounts)
    ## The remaining mutations that are found in singletons are only found in S samples and so get classified as S
    currentSetS %<>% c(setdiff(currentSet, currentSetR))
    ## Record the classification of the newly classified mutations as S or R, accordingly; leave the rest alone
    output[currentSetS, "class"] %<>% ifelse(is.na(.), SCODE, .)
    output[currentSetR, "class"] %<>% ifelse(is.na(.), RCODE, .)
    ## Record the iteration of the newly classified mutations, but leave the rest alone
    output[c(currentSetR, currentSetS), "iter"] %<>% ifelse(is.na(.), iter, .)
    ## If required, record the isolates that are SOLO for the newly classified mutations as such; leave the rest alone
    if (listIsolates) {
      isoList = isoList %>% 
        bind_rows(currentEdges %>% 
                    set_colnames(c("sample_id", "variant")) %>% 
                    as_tibble() %>% 
                    dplyr::filter(is.na(output[variant, "Tcnt"])) %>% 
                    mutate(iter = iter))
    }
    ## Record the SOLO counts of the newly classified mutations as the current SOLO counts; leave the rest alone
    output[currentSet,  "Tcnt"] %<>% ifelse(is.na(.), currentTCounts, .)
    output[currentSetR, "Rcnt"] %<>% ifelse(is.na(.), currentRCounts, .)
    ## Increment iteration counter
    iter %<>% add(1)
  }
  ## Convert the output into a tibble, replace NA values with 0 and add the original mutation names as 1st column
  initMutations = rownames(output)
  output %<>%
    as_tibble() %>%
    mutate_all(~{replace_na(., 0) %>% as.integer()}) %>%
    mutate(mutation = initMutations, .before = 1)
  output = list(class = output)
  if (listIsolates) {
     output %<>% c(list(isolate = isoList))
  }
  output
}

#' Run the full SOLO pipeline with preprocessing, maxIter steps, and postprocessing
#' @param inputTab A data frame with columns drug, sample_id, variant, phenotype, and optionally a fifth column for SOnly or het information.
#' @param maxIter Maximum number of SOLO algorithm iterations; use Inf for no limit.
#' @param stage Integer stage label appended to the output; NULL omits it.
#' @param removeSOnly Logical; if TRUE, variants occurring only in susceptible isolates are classified and removed before the main loop.
#' @inheritParams mainDriver
#' @export
runSOLOPipeline = function(inputTab, maxIter, stage = NULL, removeSOnly = TRUE, listIsolates = FALSE) {
  stopifnot(ncol(inputTab) %in% 4:5 && all(colnames(inputTab)[1:4] == c("drug", "sample_id", "variant", "phenotype")))
  ## Initialise the output table
  outputs  = tibble(drug = character(), mutation = character(), class = integer(), iter = integer(), Tcnt = integer(), Rcnt = integer())
  isoLists = tibble()
  ## For each drug, process the first maxIter steps of the original solo algorithm
  for (drugName in unique(inputTab$drug)) {
    curTab = inputTab %>%
      dplyr::filter(drug == drugName) %>%
      select(-drug) %>%
      dplyr::filter(!is.na(variant)) ## Avoids warnings about phenotype-only samples
    curOutput = originalSOLO(inputTab = curTab, maxIter = maxIter, removeSOnly = removeSOnly, listIsolates = listIsolates)
    outputs %<>%
      bind_rows(curOutput$class %>% mutate(drug = drugName, .before = 1))
    if (listIsolates) {
      isoLists %<>%
        bind_rows(curOutput$isolate %>% mutate(drug = drugName))
    }
  }
  ## Update the outputs for compatibility with the original format and return them
  outputs %<>%
    rename(variant = "mutation") %>%
    mutate_at("class", ~{recode(as.factor(.), !!!CODE_KEY) %>% as.character()}) %>%
    mutate(Scnt = Tcnt - Rcnt, Tcnt = NULL)
  if (!is.null(stage)) {
    outputs %<>%
      mutate(stage = stage)
    isoLists %<>%
      mutate(stage = stage)
  }
  outputs = list(classes = outputs)
  if (listIsolates) {
    outputs %<>% c(list(isolates = isoLists %>% mutate_at('sample_id', as.integer)))
  }
  outputs
}
