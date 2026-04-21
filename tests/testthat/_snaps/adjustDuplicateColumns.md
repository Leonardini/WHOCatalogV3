# prints warning for conflicting non-NA entries when warn = TRUE

    Code
      adjustDuplicateColumns(df, warn = TRUE)
    Output
      [1] "Warning: column col has 1 conflicting entries!"
      # A tibble: 2 x 1
          col
        <dbl>
      1     1
      2     2

