# runGenotypeConsistencyTests errors when pncA row has wrong drug

    Code
      runGenotypeConsistencyTests(df)
    Condition
      Error in `runGenotypeConsistencyTests()`:
      ! all(allGenotypes %>% dplyr::filter(gene == PZA_GENE) %>% pull(drug) ==  .... is not TRUE

# runGenotypeConsistencyTests errors when rpoB row has wrong drug

    Code
      runGenotypeConsistencyTests(df)
    Condition
      Error in `runGenotypeConsistencyTests()`:
      ! all(allGenotypes %>% dplyr::filter(gene == RRDR_GENE) %>% pull(drug) ==  .... is not TRUE

# runGenotypeConsistencyTests errors when mutation is missing but effect is not

    Code
      runGenotypeConsistencyTests(df)
    Condition
      Error in `runGenotypeConsistencyTests()`:
      ! all((allGenotypes$mutation == MISSING_VARIANT) == (allGenotypes$effect ==  .... is not TRUE

# runGenotypeConsistencyTests errors when effect is inconsistent for the same variant

    Code
      runGenotypeConsistencyTests(df)
    Condition
      Error in `runGenotypeConsistencyTests()`:
      ! all(testConsistent(allGenotypes, c("drug", "variant"), consistentVars = c("effect",  .... is not TRUE

# runGenotypeConsistencyTests errors when tier is inconsistent for the same drug-gene

    Code
      runGenotypeConsistencyTests(df)
    Condition
      Error in `runGenotypeConsistencyTests()`:
      ! all(testConsistent(allGenotypes, c("drug", "gene"), consistentVars = "tier")[[1]]) is not TRUE

