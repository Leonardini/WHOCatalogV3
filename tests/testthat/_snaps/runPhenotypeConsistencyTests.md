# runPhenotypeConsistencyTests errors when phenotype is inconsistent for same sample-category-drug

    Code
      runPhenotypeConsistencyTests(df)
    Condition
      Error in `runPhenotypeConsistencyTests()`:
      ! testConsistent(allPhenotypes, c("sample_id", "category_phenotype",  .... is not TRUE

# runPhenotypeConsistencyTests errors when same sample-drug appears in both WHO and ALL

    Code
      runPhenotypeConsistencyTests(df)
    Condition
      Error in `runPhenotypeConsistencyTests()`:
      ! anyDuplicated(allPhenotypes %>% dplyr::filter(category_phenotype %in%  .... is not TRUE

