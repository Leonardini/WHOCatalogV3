# marks variants with MAF below minMAF as het when lowMAFHet = TRUE

    Code
      applyThresholds(df, minMAF = 0.75, lowMAFHet = TRUE, minQ = NA)
    Condition
      Warning in `applyThresholds()`:
      no max(quality) column found in input data; skipping quality filtering.
    Output
      # A tibble: 2 x 3
        `max(af)` variant het  
            <dbl> <chr>   <lgl>
      1       0.9 A       FALSE
      2       0.5 B       TRUE 

# marks variant named 'missing' as het regardless of MAF

    Code
      applyThresholds(df, minMAF = 0.75, lowMAFHet = TRUE, minQ = NA)
    Condition
      Warning in `applyThresholds()`:
      no max(quality) column found in input data; skipping quality filtering.
    Output
      # A tibble: 2 x 3
        `max(af)` variant het  
            <dbl> <chr>   <lgl>
      1       0.9 A       FALSE
      2       0.9 missing TRUE 

# removes variants with MAF below minMAF when lowMAFHet = FALSE

    Code
      applyThresholds(df, minMAF = 0.75, lowMAFHet = FALSE, minQ = NA)
    Condition
      Warning in `applyThresholds()`:
      no max(quality) column found in input data; skipping quality filtering.
    Output
      # A tibble: 1 x 2
        `max(af)` variant
            <dbl> <chr>  
      1       0.9 A      

# removes variants with quality below minQ when lowQHet = FALSE

    Code
      applyThresholds(df, minMAF = NA, minQ = 1000, lowQHet = FALSE)
    Condition
      Warning in `applyThresholds()`:
      no max(af) column found in input data; skipping MAF filtering.
    Output
      # A tibble: 1 x 3
        `max(af)` `max(quality)` variant
            <dbl>          <dbl> <chr>  
      1       0.9           1200 A      

# skips MAF filtering and prints warning when max(af) column is absent

    Code
      applyThresholds(df, minMAF = 0.75, minQ = NA)
    Condition
      Warning in `applyThresholds()`:
      no max(af) column found in input data; skipping MAF filtering.
      Warning in `applyThresholds()`:
      no max(quality) column found in input data; skipping quality filtering.
    Output
      # A tibble: 2 x 1
        variant
        <chr>  
      1 A      
      2 B      

# skips quality filtering and prints warning when max(quality) column is absent

    Code
      applyThresholds(df, minMAF = NA, minQ = 1000)
    Condition
      Warning in `applyThresholds()`:
      no max(af) column found in input data; skipping MAF filtering.
      Warning in `applyThresholds()`:
      no max(quality) column found in input data; skipping quality filtering.
    Output
      # A tibble: 2 x 2
        `max(af)` variant
            <dbl> <chr>  
      1       0.9 A      
      2       0.9 B      

