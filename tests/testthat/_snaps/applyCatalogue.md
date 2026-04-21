# applyCatalogue errors when required columns are missing

    Code
      applyCatalogue(tibble(drug = "Isoniazid"), tf, minMAF = NA, minQ = NA)
    Condition
      Error in `applyCatalogue()`:
      ! Input data must contain 'drug', 'effect', 'gene', 'mutation' and 'variant' columns.

