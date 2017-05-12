# methylation_EWAS_workshop

This is a basic EWAS tutorial, introducing concepts for analysing methylation data.

## Data

The workshop uses public Illumina 450K data from several studies (GSE40279, GSE48472, GSE41114, GSE42700, GSE46573, GSE48472, GSE50586). It is currently designed to investigate the methylation differences between blood (whole blood) and buccal (cheek swabs) tissues. It is aimed at giving a very basic overview of the structure and type of data used in EWAS as well as a few of the available tools and methods to explore this data.

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

You can download this repository as a zip file (see top of this page for Download button) and extract to the folder `methylation_EWAS_workshop` where appropriate on your system.

**Clone**

If you are familiar with git and GitHub you can clone into this repository:

`git clone https://github.com/sirselim/methylation_EWAS_workshop.git`

### 2. move to the correct directory and start working

Open `R`/`RStudio` and ensure that your working directory is set to the `methylation_EWAS_workshop` directory, then open the `EWAS_analysis.R` script and follow along.

