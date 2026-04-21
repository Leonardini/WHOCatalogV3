# errors when manualTab has fewer than 2 columns

    Code
      recodeValues(tibble(key = "a", val = 1), tibble(key = "a"))
    Condition
      Error in `recodeValues()`:
      ! ncol(manualTab) == 2 is not TRUE

# errors when manualTab has more than 2 columns

    Code
      recodeValues(tibble(key = "a", val = 1), tibble(a = 1, b = 2, c = 3))
    Condition
      Error in `recodeValues()`:
      ! ncol(manualTab) == 2 is not TRUE

# errors when a manualTab column name ends with _first or _second

    Code
      recodeValues(tibble(key = "a", val = 1), tibble(key = "a", val_first = 10))
    Condition
      Error in `recodeValues()`:
      ! all(str_ends(coln, "_first", negate = TRUE)) && all(str_ends(coln,  .... is not TRUE

