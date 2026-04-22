# adjustStatsForRemovedMutations prints a message

    Code
      adjustStatsForRemovedMutations(stats, dataset, remove, 2L, make_empty_solos())
    Message
      Removing 4 isolates at stage 2
    Output
      # A tibble: 1 x 21
        drug   variant absent_R absent_S SOLO_SorR stage  tier neutral datasets effect
        <chr>  <chr>      <int>    <int>     <dbl> <int> <int> <lgl>   <chr>    <chr> 
      1 Isoni~ katG_p~        4        5         2     2     1 FALSE   ALL      misse~
      # i 11 more variables: position <chr>, lit_mutation <lgl>, prev_version <lgl>,
      #   correctAll <lgl>, correctSOLO <lgl>, pos1 <int>, present <int>,
      #   present_R <int>, present_S <int>, SOLO_R <dbl>, SOLO_S <dbl>

