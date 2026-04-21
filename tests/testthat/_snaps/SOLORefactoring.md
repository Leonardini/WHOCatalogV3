# runSOLOPipeline errors when required column names are wrong

    Code
      runSOLOPipeline(tibble(a = 1, b = 2, c = 3, d = 4), maxIter = 1L)
    Condition
      Error in `runSOLOPipeline()`:
      ! ncol(inputTab) %in% 4:5 && all(colnames(inputTab)[1:4] == c("drug",  .... is not TRUE

