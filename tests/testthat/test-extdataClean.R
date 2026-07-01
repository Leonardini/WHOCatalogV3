testthat::local_edition(3)

read_text_file <- function(f) {
  raw <- readBin(f, "raw", n = file.info(f)$size)
  s <- rawToChar(raw)
  if (validUTF8(s)) Encoding(s) <- "UTF-8" else s <- iconv(s, from = "latin1", to = "UTF-8")
  s
}

JOIN_BREAKING_CODEPOINTS <- c(
  setdiff(0x00:0x1F, c(0x0A, 0x0D)),
  0x7F, 0x80:0x9F,
  0xA0, 0xAD, 0x061C, 0x1680, 0x180E,
  0x2000:0x200F, 0x2028, 0x2029, 0x202A:0x202F,
  0x205F, 0x2060, 0x2066:0x2069, 0x3000, 0xFEFF, 0xFFF9:0xFFFB
)

test_that("bundled extdata CSVs contain no hidden characters that could break joins or merges", {
  extdata <- system.file("extdata", package = "SOLOport")
  skip_if(extdata == "", "extdata directory not found")
  csvFiles <- list.files(extdata, pattern = "\\.csv$", full.names = TRUE)
  skip_if(length(csvFiles) == 0, "no extdata CSV files found")

  offenders <- lapply(csvFiles, function(f) {
    bad <- sort(unique(intersect(utf8ToInt(read_text_file(f)), JOIN_BREAKING_CODEPOINTS)))
    if (length(bad)) sprintf("%s: %s", basename(f), paste(sprintf("U+%04X", bad), collapse = ", ")) else NULL
  })
  expect_equal(Filter(Negate(is.null), offenders), list())
})
