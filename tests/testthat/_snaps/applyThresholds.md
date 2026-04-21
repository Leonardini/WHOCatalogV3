# skips MAF filtering and prints warning when max(af) column is absent

    Code
      applyThresholds(df, minMAF = 0.75, minQ = NA)
    Output
      [1] "Warning: no max(af) column found in input data; skipping MAF filtering."
      [1] "Warning: no max(quality) column found in input data; skipping quality filtering."
      # A tibble: 2 x 1
        variant
        <chr>  
      1 A      
      2 B      

# skips quality filtering and prints warning when max(quality) column is absent

    Code
      applyThresholds(df, minMAF = NA, minQ = 1000)
    Output
      [1] "Warning: no max(af) column found in input data; skipping MAF filtering."
      [1] "Warning: no max(quality) column found in input data; skipping quality filtering."
      # A tibble: 2 x 2
        `max(af)` variant
            <dbl> <chr>  
      1       0.9 A      
      2       0.9 B      

