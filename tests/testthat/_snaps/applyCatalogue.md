# applyCatalogue errors when required columns are missing

    Code
      applyCatalogue(tibble(drug = "Isoniazid"), tf, minMAF = NA, minQ = NA)
    Condition
      Error in `applyCatalogue()`:
      ! Input data must contain 'drug', 'effect', 'gene', 'mutation' and 'variant' columns.

# applyCatalogue assigns the catalogue grade to a matching variant

    Code
      applyCatalogue(make_ac_row("Isoniazid", "katG", "p.Ser315Thr"), tf, minMAF = NA,
      minQ = NA, LoF = FALSE)
    Condition
      Warning in `applyThresholds()`:
      no max(af) column found in input data; skipping MAF filtering.
      Warning in `applyThresholds()`:
      no max(quality) column found in input data; skipping quality filtering.
    Output
      # A tibble: 1 x 9
        drug      gene  variant      mutation effect  pos1  pos2 RRDR_NON_SILENT Final
        <chr>     <chr> <chr>        <chr>    <chr>  <int> <int> <lgl>           <int>
      1 Isoniazid katG  katG_p.Ser3~ p.Ser31~ misse~   315    NA FALSE               1

# applyCatalogue leaves Final as NA for variants not in the catalogue

    Code
      applyCatalogue(make_ac_row("Isoniazid", "katG", "p.Arg463Leu"), tf, minMAF = NA,
      minQ = NA, LoF = FALSE)
    Condition
      Warning in `applyThresholds()`:
      no max(af) column found in input data; skipping MAF filtering.
      Warning in `applyThresholds()`:
      no max(quality) column found in input data; skipping quality filtering.
    Output
      # A tibble: 1 x 9
        drug      gene  variant      mutation effect  pos1  pos2 RRDR_NON_SILENT Final
        <chr>     <chr> <chr>        <chr>    <chr>  <int> <int> <lgl>           <int>
      1 Isoniazid katG  katG_p.Arg4~ p.Arg46~ misse~   463    NA FALSE              NA

# applyCatalogue sets Final = 2 for unmatched RRDR non-silent variants

    Code
      applyCatalogue(make_ac_row("Rifampicin", "rpoB", "p.Ser450Leu"), tf, minMAF = NA,
      minQ = NA, LoF = FALSE)
    Condition
      Warning in `applyThresholds()`:
      no max(af) column found in input data; skipping MAF filtering.
      Warning in `applyThresholds()`:
      no max(quality) column found in input data; skipping quality filtering.
    Output
      # A tibble: 1 x 9
        drug       gene  variant     mutation effect  pos1  pos2 RRDR_NON_SILENT Final
        <chr>      <chr> <chr>       <chr>    <chr>  <int> <int> <lgl>           <int>
      1 Rifampicin rpoB  rpoB_p.Ser~ p.Ser45~ misse~   450    NA TRUE                2

# applyCatalogue sets Final = 2 for unmatched LoF variants when pooled LoF is graded 1 or 2

    Code
      applyCatalogue(make_ac_row("Isoniazid", "katG", "p.Ala12Val", effect = "frameshift"),
      tf, minMAF = NA, minQ = NA, LoF = TRUE)
    Condition
      Warning in `applyThresholds()`:
      no max(af) column found in input data; skipping MAF filtering.
      Warning in `applyThresholds()`:
      no max(quality) column found in input data; skipping quality filtering.
    Output
      # A tibble: 1 x 10
        drug      gene  variant      mutation effect  pos1  pos2 RRDR_NON_SILENT Final
        <chr>     <chr> <chr>        <chr>    <chr>  <int> <int> <lgl>           <int>
      1 Isoniazid katG  katG_p.Ala1~ p.Ala12~ frame~    12    NA FALSE               2
      # i 1 more variable: LoF_candidate <lgl>

# applyCatalogue does not apply the LoF rule when LoF = FALSE

    Code
      applyCatalogue(make_ac_row("Isoniazid", "katG", "p.Ala12Val", effect = "frameshift"),
      tf, minMAF = NA, minQ = NA, LoF = FALSE)
    Condition
      Warning in `applyThresholds()`:
      no max(af) column found in input data; skipping MAF filtering.
      Warning in `applyThresholds()`:
      no max(quality) column found in input data; skipping quality filtering.
    Output
      # A tibble: 1 x 9
        drug      gene  variant      mutation effect  pos1  pos2 RRDR_NON_SILENT Final
        <chr>     <chr> <chr>        <chr>    <chr>  <int> <int> <lgl>           <int>
      1 Isoniazid katG  katG_p.Ala1~ p.Ala12~ frame~    12    NA FALSE              NA

