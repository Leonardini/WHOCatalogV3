source("SOLOport/Driver.R")

last_extract = "2023-04-25T06_00_10.443990_jr_b741dc136e079fa8583604a4915c0dc751724ae9880f06e7c2eacc939e086536"

mainDriver(
    NON_DATABASE_DIRECTORY = Sys.getenv("CODEBUILD_SRC_DIR_DTA"),
    DATA_DIRECTORY = Sys.getenv("CODEBUILD_SRC_DIR_S3"),
    EXTRACTION_ID = last_extract,
    OUTPUT_DIRECTORY = paste0("Results/", last_extract)
)