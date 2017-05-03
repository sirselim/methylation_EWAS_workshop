# methylation_EWAS_workshop
This is a basic EWAS tutorial, introducing basic concepts for analysing methylation data.

## Data
The workshop uses public Illumina 450K data from several studies (GSE40279, GSE48472), and is currently designed to investigate the methylation differences between blood (whole blood) and buccal (cheek swabs) tissues.

## Required software and packages

You will need to ensure you have the following installed

  - the latest version of `R`: https://cran.r-project.org/
  - `RStudio` (optional but recommended): https://www.rstudio.com/products/rstudio/download2/
  - the `R` package `minfi`: https://www.bioconductor.org/packages/release/bioc/html/minfi.html
    + to install `minfi` run the following within `R`/`RStudio`:
      - `source("https://bioconductor.org/biocLite.R")`
      - `biocLite("minfi")`

*NOTE: this workshop should work across all OS (Windows, Linux, Mac).*
