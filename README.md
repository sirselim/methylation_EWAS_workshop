# methylation_EWAS_workshop

This is a basic EWAS tutorial, introducing basic concepts for analysing methylation data.

## Data

The workshop uses public Illumina 450K data from several studies (GSE40279, GSE48472), and is currently designed to investigate the methylation differences between blood (whole blood) and buccal (cheek swabs) tissues.

## Required software and packages

You will need to ensure you have the following installed on your system:

  - the latest version of `R`: https://cran.r-project.org/
  - `RStudio` (optional but recommended): https://www.rstudio.com/products/rstudio/download2/
  - the `R` package `minfi`: https://www.bioconductor.org/packages/release/bioc/html/minfi.html
    + to install `minfi` run the following within `R`/`RStudio`:
      - `source("https://bioconductor.org/biocLite.R")`
      - `biocLite("minfi")`

*NOTE: this workshop should work across all OS (Windows, Linux, Mac) as long as the above are correctly installed.*

## Set up

### 1. obtain the workshop files

Once the above software and R packages are installed set up is as simple as downloading or cloning this repository, it contains all the data and scripts required for the workshop.

**Download**

You can download this repository as a zip file and extract to the folder `methylation_EWAS_workshop` where appropriate on your system.

**Clone**

If you are familiar with git and GitHub you can clone into this repository:

`git clone https://github.com/sirselim/methylation_EWAS_workshop.git`

### 2. move to the correct directory and start working

Open `R`/`RStudio` and ensure that your working directory is set to the `methylation_EWAS_workshop` directory, then open the `EWAS_analysis.R` script and follow along.

