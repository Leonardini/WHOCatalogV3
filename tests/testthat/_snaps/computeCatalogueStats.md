# computeCatalogueStats prints the number of rows being processed

    Code
      computeCatalogueStats(make_stats_df())
    Output
      [1] "There are 1 rows to process"
      # A tibble: 1 x 51
        drug  variant present_R present_S absent_R absent_S SOLO_R SOLO_S correctAll
        <chr> <chr>       <int>     <int>    <int>    <int>  <int>  <int> <lgl>     
      1 DrugA V1             10         2        3       15      8      1 TRUE      
      # i 42 more variables: correctSOLO <lgl>, PPV <dbl>, PPV_SOLO <dbl>,
      #   PPVc_SOLO <dbl>, Sens <dbl>, Sens_SOLO <dbl>, Spec <dbl>, Spec_SOLO <dbl>,
      #   PPV_lb <dbl>, PPV_ub <dbl>, PPV_SOLO_lb <dbl>, PPV_SOLO_ub <dbl>,
      #   PPVc_SOLO_lb <dbl>, PPVc_SOLO_ub <dbl>, Sens_lb <dbl>, Sens_ub <dbl>,
      #   Sens_SOLO_lb <dbl>, Sens_SOLO_ub <dbl>, Spec_lb <dbl>, Spec_ub <dbl>,
      #   Spec_SOLO_lb <dbl>, Spec_SOLO_ub <dbl>, OR <dbl>, OR_SOLO <dbl>,
      #   OR_exact_lb <dbl>, OR_exact_ub <dbl>, OR_pvalue <dbl>, ...

