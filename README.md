# SOLOport

An R package implementing the SOLO algorithm for analysing genomic variants and their association with drug resistance phenotypes, used to produce the WHO tuberculosis mutation catalogue.

## Installation

```r
# install.packages("remotes")
remotes::install_github("lchindel/SOLOport")
```

## Dependencies

```r
install.packages(c("dplyr", "igraph", "magrittr", "purrr", "readr", "readxl",
                   "stringr", "tibble", "tidyr"))
```

## Usage

### Full analysis pipeline

```r
library(SOLOport)

mainDriver(
  EXTRACTION_ID = "<your-extraction-id>",
  OUTPUT_DIRECTORY = "Results/<your-extraction-id>",
  DATA_DIRECTORY  = "path/to/database/extraction/files/"
)
```

`mainDriver()` orchestrates the complete analysis: genotype preprocessing, sample QC, the iterative SOLO pipeline, neutral variant identification, mutation grading, and sensitivity/specificity computation.

### Core functions

| Function | Purpose |
|---|---|
| `runSOLOPipeline()` | Run the iterative SOLO algorithm on a genotype-phenotype table |
| `neutralAlgorithm()` | Identify neutral variants using a PPV-based multi-set algorithm |
| `gradeMutations()` | Apply grading rules and write the graded catalogue |
| `applyGradingRules()` | Pure-computation grading (no file I/O) |
| `computeSensSpec()` | Compute sensitivity/specificity statistics for a catalogue |
| `applyCatalogue()` | Apply a WHO catalogue to a genotype-phenotype dataset |
| `computeCatalogueStats()` | Compute per-drug catalogue statistics |
| `computeDominantLineage()` | Assign dominant lineage to each sample |
| `mainDriver()` | End-to-end orchestrator |

### Bundled reference data

The package ships with the following reference files in `inst/extdata/`:

- `who_catalogue_prev_version.xlsx` — previous-version WHO catalogue used by the neutral algorithm
- `literature_neutrals.csv` — literature-derived neutral mutations
- `assay_mutations.csv` — assay-related mutations
- `allelic_exchanges.csv` — allelic exchange experiment results
- `grading_comments.csv` — grading comment categories
- `manual_check.csv` — manually reviewed variant calls

## Running tests

```r
devtools::load_all()
testthat::test_dir("tests/testthat")
```

## License

MIT © Leonid Chindelevitch. See [LICENSE.md](LICENSE.md).
